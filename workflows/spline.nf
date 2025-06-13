#!/usr/bin/env nextflow

// Import modules
include { TRENDY } from '../modules/local/trendy'
include { PREPROCESS } from '../modules/local/preprocess'
include { ANALYZE_DATA_SIZE } from '../modules/local/analyze_data_size'

workflow SPLINE {
    // Define input channels
    if (params.pathogen) {
        // Convert comma-separated string to a channel of pathogens
        def pathogenList = params.pathogen.tokenize(',')
        pathogens = Channel.fromList(pathogenList)
    } else {
        // Default to CAMPYLOBACTER and CYCLOSPORA for testing
        pathogens = Channel.of('CAMPYLOBACTER', 'CYCLOSPORA')
    }
    
    // Handle pathogen grouping if specified
    if (params.pathogen_grouping) {
        // Parse pathogen:subgroup format using pipe delimiter
        def groupingList = params.pathogen_grouping.tokenize('|')
        pathogenGrouping = Channel.fromList(groupingList)
    } else {
        // If no grouping specified, use simple pathogen list
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
        dataMetrics = ANALYZE_DATA_SIZE.out.metrics
            .map { metrics_file ->
                def json_text = metrics_file.text
                def metrics = new groovy.json.JsonSlurper().parseText(json_text)
                return metrics
            }

        // Parse pathogen groupings and combine with metrics
        pathogenGroupingWithMetrics = pathogenGrouping
            .map { grouping ->
                // Parse pathogen:subgroup format
                def parts = grouping.split(':')
                def pathogen = parts[0]
                def subgroup = parts.length > 1 ? parts[1] : 'combined'
                tuple(grouping, pathogen, subgroup)
            }
            .combine(dataMetrics)
            .map { grouping, pathogen, subgroup, metrics ->
                tuple(grouping, pathogen, subgroup, metrics[pathogen] ?: [rows: 0, complexity: 0])
            }

        // Run TRENDY with preprocessed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics.map { it[0] }, // full grouping string (e.g., "SALMONELLA:Enteritidis")
            pathogenGroupingWithMetrics.map { it[1] }, // base pathogen (e.g., "SALMONELLA")
            pathogenGroupingWithMetrics.map { it[2] }, // subgroup (e.g., "Enteritidis")
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.projID,
            params.trendyScript,
            params.preprocessed,
            cleanFile,
            pathogenGroupingWithMetrics.map { it[3] },  // metrics
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
        dataMetrics = ANALYZE_DATA_SIZE.out.metrics
            .map { metrics_file ->
                def json_text = metrics_file.text
                def metrics = new groovy.json.JsonSlurper().parseText(json_text)
                return metrics
            }

        // Parse pathogen groupings and combine with metrics
        pathogenGroupingWithMetrics = pathogenGrouping
            .map { grouping ->
                // Parse pathogen:subgroup format
                def parts = grouping.split(':')
                def pathogen = parts[0]
                def subgroup = parts.length > 1 ? parts[1] : 'combined'
                tuple(grouping, pathogen, subgroup)
            }
            .combine(dataMetrics)
            .map { grouping, pathogen, subgroup, metrics ->
                tuple(grouping, pathogen, subgroup, metrics[pathogen] ?: [rows: 0, complexity: 0])
            }

        // Run TRENDY with processed data and metrics
        TRENDY(
            pathogenGroupingWithMetrics.map { it[0] }, // full grouping string
            pathogenGroupingWithMetrics.map { it[1] }, // base pathogen
            pathogenGroupingWithMetrics.map { it[2] }, // subgroup
            mmwrFile,
            censusFileB,
            censusFileP,
            params.travel,
            params.cidt,
            params.projID,
            params.trendyScript,
            true,
            processedFile,
            pathogenGroupingWithMetrics.map { it[3] },  // metrics
            catchmentConfig
        )
    }

    // Log completion
    log.info "Pipeline completed successfully"
}
