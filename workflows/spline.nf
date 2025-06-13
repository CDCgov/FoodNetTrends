#!/usr/bin/env nextflow

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'
include { ANALYZE_DATA_SIZE } from '../modules/local/analyze_data_size'

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

        // Analyze data size for resource allocation
        ANALYZE_DATA_SIZE(cleanFile)
        
        // Read the metrics and create a map
        metricsChannel = ANALYZE_DATA_SIZE.out.metrics
            .map { metrics_file ->
                def json_text = metrics_file.text
                def metrics = new groovy.json.JsonSlurper().parseText(json_text)
                return metrics
            }
        
        // Parse pathogen groupings and combine with metrics
        // Handle case where pathogenGrouping might be null (AUTO_DISCOVER)
        if (!pathogenGrouping) {
            error "Pathogen grouping is not defined. This should not happen for preprocessed data."
        }
        
        pathogenGroupingWithMetrics = pathogenGrouping
            .combine(metricsChannel)
            .map { grouping, metrics ->
                // Parse pathogen:subgroup format
                def parts = grouping.split(':')
                def pathogen = parts[0]
                def subgroup = parts.length > 1 ? parts[1] : 'combined'
                def pathogenMetrics = metrics[pathogen] ?: [rows: 0, complexity: 0]
                if (pathogenMetrics.rows == 0) {
                    log.warn "No data found for pathogen: ${pathogen}. Using default metrics."
                }
                tuple(grouping, pathogen, subgroup, pathogenMetrics)
            }

        // Run TRENDY with preprocessed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
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

        // Analyze data size for resource allocation
        ANALYZE_DATA_SIZE(processedFile)
        
        // Read the metrics and create a map
        metricsChannel = ANALYZE_DATA_SIZE.out.metrics
            .map { metrics_file ->
                def json_text = metrics_file.text
                def metrics = new groovy.json.JsonSlurper().parseText(json_text)
                return metrics
            }
        
        // If AUTO_DISCOVER, extract pathogens from the metrics
        if (params.pathogen == 'AUTO_DISCOVER') {
            // Extract pathogen list from metrics and create groupings
            pathogenGrouping = metricsChannel
                .map { metrics ->
                    def pathogenList = metrics.keySet().toList()
                    if (pathogenList.isEmpty()) {
                        error "No pathogens found in preprocessed data. Check if preprocessing completed successfully."
                    }
                    log.info "Auto-discovered pathogens: ${pathogenList.join(', ')}"
                    return pathogenList
                }
                .flatMap()
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
                def pathogenMetrics = metrics[pathogen] ?: [rows: 0, complexity: 0]
                if (pathogenMetrics.rows == 0) {
                    log.warn "No data found for pathogen: ${pathogen}. Using default metrics."
                }
                tuple(grouping, pathogen, subgroup, pathogenMetrics)
            }

        // Run TRENDY with processed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics,
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
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
