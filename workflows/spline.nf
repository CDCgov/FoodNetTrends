#!/usr/bin/env nextflow

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'
include { RESOURCE_PROFILER } from '../modules/local/resource_profiler'

workflow SPLINE {
    // Define input channels
    if (params.pathogen && params.pathogen != 'AUTO_DISCOVER') {
        // Convert comma-separated string to a channel of pathogens
        def pathogenList = params.pathogen.tokenize(',')
        pathogens = Channel.fromList(pathogenList)
    } else if (params.pathogen == 'AUTO_DISCOVER') {
        // Will be populated after preprocessing
        pathogens = null
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
    }
    
    // Handle pathogen grouping if specified
    if (params.pathogen_grouping && params.pathogen_grouping.trim()) {
        // Parse pathogen:subgroup format using pipe delimiter
        def groupingList = params.pathogen_grouping.tokenize('|')
        pathogenGrouping = Channel.fromList(groupingList)
    } else if (params.pathogen == 'AUTO_DISCOVER') {
        // Defer pathogenGrouping creation until after preprocessing
        pathogenGrouping = null
    } else if (pathogens) {
        // If no grouping specified, use simple pathogen list
        pathogenGrouping = pathogens.map { p -> "${p}:combined" }
    } else {
        // Fallback to default pathogens if nothing specified
        log.warn "No pathogens specified, using defaults: CAMPYLOBACTER, CYCLOSPORA"
        pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
        pathogenGrouping = pathogens.map { p -> "${p}:combined" }
    }

    // Input files
    mmwrFile = file(params.mmwrFile)
    censusFileB = file(params.censusFileB)
    censusFileP = file(params.censusFileP)
    
    // Configuration files (optional)
    serotypeConfig = params.serotype_config ? file(params.serotype_config) : file('NO_FILE')
    catchmentConfig = params.catchment_config ? file(params.catchment_config) : file('NO_FILE')

    // Check if files exist
    if (!mmwrFile.exists()) {
        error "MMWR file not found: ${params.mmwrFile}"
    }
    if (!censusFileB.exists()) {
        error "Census bacterial file not found: ${params.censusFileB}"
    }
    if (!censusFileP.exists()) {
        error "Census parasitic file not found: ${params.censusFileP}"
    }

    // Log pipeline start
    log.info """
    ==============================================
    FoodNet Trends Pipeline
    ==============================================
    Project ID    : ${params.projID}
    MMWR File     : ${params.mmwrFile}
    Census Files  : ${params.censusFileB}, ${params.censusFileP}
    States        : ${params.states ?: 'all'}
    Travel        : ${params.travel}
    CIDT          : ${params.cidt}
    Pathogens     : ${params.pathogen ?: 'default (CAMPYLOBACTER,CYCLOSPORA)'}
    Cores         : ${params.cpus ?: 'default'}
    Chains        : ${params.chains}
    Iterations    : ${params.iterations}
    Adapt Delta   : ${params.adapt_delta}
    Max Treedepth : ${params.max_treedepth}
    Seed          : ${params.seed}
    Output Dir    : ${params.outdir}/${params.projID}
    ==============================================
    """

    // Conditional preprocessing
    if (params.preprocessed) {
        log.info "Using preprocessed data from: ${params.cleanFile}"
        cleanFile = file(params.cleanFile)
        if (!cleanFile.exists()) {
            error "Preprocessed file not found: ${params.cleanFile}"
        }

        // Check for existing resource profile
        def resourceProfilePath = cleanFile.parent.resolve("resource_profile.csv")
        def resourceProfile = file(resourceProfilePath)
        
        if (resourceProfile.exists()) {
            log.info "Using existing resource profile: ${resourceProfilePath}"
            // Read the CSV directly
            metricsChannel = Channel.fromPath(resourceProfilePath)
                .splitCsv(header: true, sep: ',', strip: true)
                .map { row ->
                    // Each row is a map with the CSV columns
                    tuple(row.pathogen, [
                        rows: row.rows as Integer,
                        sites: row.sites as Integer,
                        years: row.years as Integer,
                        complexity: row.complexity as Long,
                        size_category: row.size_category
                    ])
                }
                .collect() // Collect all tuples into a list
                .map { tupleList ->
                    // Convert list of tuples into a map
                    if (tupleList.isEmpty()) {
                        error "Resource profile is empty - no pathogen data found"
                    }
                    def metrics = [:]
                    tupleList.each { tuple ->
                        metrics[tuple[0]] = tuple[1]
                    }
                    return metrics
                }
        } else {
            log.info "Generating resource profile for preprocessed data"
            // Run resource profiler
            RESOURCE_PROFILER(cleanFile)
            
            // Read the CSV output
            metricsChannel = RESOURCE_PROFILER.out.profile
                .map { csvFile -> 
                    // Read the CSV file content and parse it
                    def metrics = [:]
                    csvFile.splitCsv(header: true, sep: ',', strip: true).each { row ->
                        log.debug "CSV row: ${row}"
                        // Convert row values to appropriate types
                        metrics[row.pathogen] = [
                            rows: row.rows as Integer,
                            sites: row.sites as Integer, 
                            years: row.years as Integer,
                            complexity: row.complexity as Long,
                            size_category: row.size_category
                        ]
                    }
                    if (metrics.isEmpty()) {
                        error "Resource profile is empty - no pathogen data found"
                    }
                    log.info "Parsed metrics for ${metrics.size()} pathogens: ${metrics.keySet().join(', ')}"
                    return metrics
                }
        }
        
        // Parse pathogen groupings and combine with metrics
        // Handle case where pathogenGrouping might be null (AUTO_DISCOVER)
        if (!pathogenGrouping) {
            // This happens when using AUTO_DISCOVER with preprocessed data
            if (params.pathogen == 'AUTO_DISCOVER') {
                log.info "Creating pathogen groupings from discovered pathogens"
                pathogenGrouping = metricsChannel
                    .flatMap { metrics ->
                        def pathogenList = metrics.keySet().toList()
                        if (pathogenList.isEmpty()) {
                            error "No pathogens found in metrics. Check if preprocessing completed successfully."
                        }
                        log.info "Auto-discovered pathogens: ${pathogenList.join(', ')}"
                        return pathogenList
                    }
                    .map { p -> "${p}:combined" }
            } else {
                error "Pathogen grouping is not defined. This should not happen for preprocessed data."
            }
        }
        
        pathogenGroupingWithMetrics = pathogenGrouping
            .combine(metricsChannel)
            .map { grouping, metrics ->
                // Parse pathogen:subgroup format
                def parts = grouping.split(':')
                def pathogen = parts[0]
                def subgroup = parts.length > 1 ? parts[1] : 'combined'
                
                // Debug: log the metrics map
                log.debug "Metrics map keys: ${metrics.keySet()}"
                log.debug "Looking for pathogen: '${pathogen}'"
                
                // Ensure we get a proper map, not just a value
                def rawMetrics = metrics[pathogen]
                def pathogenMetrics
                if (rawMetrics instanceof Map) {
                    pathogenMetrics = rawMetrics
                } else if (rawMetrics instanceof List && rawMetrics.size() > 0) {
                    // If it's a list, try to extract values
                    log.warn "Metrics for ${pathogen} is a list, not a map: ${rawMetrics}"
                    pathogenMetrics = [rows: rawMetrics[0], complexity: 0]
                } else if (rawMetrics) {
                    // Single value, assume it's rows
                    log.warn "Metrics for ${pathogen} is a single value: ${rawMetrics}"
                    pathogenMetrics = [rows: rawMetrics, complexity: 0]
                } else {
                    // No data
                    pathogenMetrics = [rows: 0, complexity: 0]
                }
                if (pathogenMetrics.rows == 0) {
                    log.warn "No data found for pathogen: ${pathogen}. Using default metrics."
                }
                // Debug: log what we're passing
                log.debug "Creating tuple for ${pathogen}: grouping=${grouping}, subgroup=${subgroup}, metrics=${pathogenMetrics}"
                tuple(grouping, pathogen, subgroup, pathogenMetrics)
            }
            .filter { grouping, pathogen, subgroup, pathogenMetrics ->
                if (pathogenMetrics.rows == 0) {
                    log.warn "Skipping ${pathogen} - no data available in preprocessed file"
                    return false
                }
                return true
            }

        // Run TRENDY with preprocessed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.states,
            params.projID,
            params.trendyScript,
            params.preprocessed,
            cleanFile,
            catchmentConfig
        )
    } else {
        log.info "Preprocessing raw data files"

        // Run preprocessing step
        PREPROCESS(
            mmwrFile,
            params.projID,
            serotypeConfig
        )

        // Create a proper channel from the preprocessed file
        processedFile = PREPROCESS.out.cleanFile

        // Generate resource profile for new data
        RESOURCE_PROFILER(processedFile)
        
        // Read the CSV output - need to read file content first
        metricsChannel = RESOURCE_PROFILER.out.profile
            .map { csvFile -> 
                // Read the CSV file content and parse it
                def metrics = [:]
                csvFile.splitCsv(header: true, sep: ',', strip: true).each { row ->
                    log.debug "CSV row: ${row}"
                    // Convert row values to appropriate types
                    metrics[row.pathogen] = [
                        rows: row.rows as Integer,
                        sites: row.sites as Integer, 
                        years: row.years as Integer,
                        complexity: row.complexity as Long,
                        size_category: row.size_category
                    ]
                }
                log.info "Parsed metrics for ${metrics.size()} pathogens: ${metrics.keySet().join(', ')}"
                return metrics
            }
        
        // If AUTO_DISCOVER, extract pathogens from the metrics
        if (params.pathogen == 'AUTO_DISCOVER') {
            // Extract pathogen list from metrics and create groupings
            pathogenGrouping = metricsChannel
                .flatMap { metrics ->
                    def pathogenList = metrics.keySet().toList()
                    if (pathogenList.isEmpty()) {
                        error "No pathogens found in preprocessed data. Check if preprocessing completed successfully."
                    }
                    log.info "Auto-discovered pathogens: ${pathogenList.join(', ')}"
                    return pathogenList
                }
                .map { p -> "${p}:combined" }
            
            // Create pathogens channel for consistency
            pathogens = pathogenGrouping.map { grouping ->
                grouping.split(':')[0]
            }
        }
        
        if (!pathogenGrouping) {
            // This handles the edge case where pathogenGrouping wasn't set earlier
            log.warn "Pathogen grouping was not properly initialized. Using defaults."
            if (!pathogens) {
                pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
            }
            pathogenGrouping = pathogens.map { p -> "${p}:combined" }
        }

        // Parse pathogen groupings and combine with metrics
        // Ensure pathogenGrouping exists before using it
        if (!pathogenGrouping) {
            error "Pathogen grouping is not defined after preprocessing. This indicates a logic error."
        }
        
        pathogenGroupingWithMetrics = pathogenGrouping
            .combine(metricsChannel)
            .map { grouping, metrics ->
                // Parse pathogen:subgroup format
                def parts = grouping.split(':')
                def pathogen = parts[0]
                def subgroup = parts.length > 1 ? parts[1] : 'combined'
                
                // Debug: log the metrics map
                log.debug "Metrics map keys: ${metrics.keySet()}"
                log.debug "Looking for pathogen: '${pathogen}'"
                
                // Ensure we get a proper map, not just a value
                def rawMetrics = metrics[pathogen]
                def pathogenMetrics
                if (rawMetrics instanceof Map) {
                    pathogenMetrics = rawMetrics
                } else if (rawMetrics instanceof List && rawMetrics.size() > 0) {
                    // If it's a list, try to extract values
                    log.warn "Metrics for ${pathogen} is a list, not a map: ${rawMetrics}"
                    pathogenMetrics = [rows: rawMetrics[0], complexity: 0]
                } else if (rawMetrics) {
                    // Single value, assume it's rows
                    log.warn "Metrics for ${pathogen} is a single value: ${rawMetrics}"
                    pathogenMetrics = [rows: rawMetrics, complexity: 0]
                } else {
                    // No data
                    pathogenMetrics = [rows: 0, complexity: 0]
                }
                if (pathogenMetrics.rows == 0) {
                    log.warn "No data found for pathogen: ${pathogen}. Using default metrics."
                }
                // Debug: log what we're passing
                log.debug "Creating tuple for ${pathogen}: grouping=${grouping}, subgroup=${subgroup}, metrics=${pathogenMetrics}"
                tuple(grouping, pathogen, subgroup, pathogenMetrics)
            }
            .filter { grouping, pathogen, subgroup, pathogenMetrics ->
                if (pathogenMetrics.rows == 0) {
                    log.warn "Skipping ${pathogen} - no data available in preprocessed file"
                    return false
                }
                return true
            }

        // Run TRENDY with processed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.states,
            params.projID,
            params.trendyScript,
            true,
            processedFile,
            catchmentConfig
        )
    }

    // Log completion
    log.info "Pipeline completed successfully"
}
