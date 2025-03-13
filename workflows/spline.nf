nextflow.enable.dsl = 2

// Create a channel with the list of pathogens to model.
// This has been modified to only run 2 pathogens for testing. It needs to have the others added before production.
Channel
    .from(['CAMPYLOBACTER', 'CYCLOSPORA'])
    .set { pathogens_ch }

// Include the modules
include { TRENDY } from '../modules/local/trendy.nf'
include { PREPROCESS } from '../modules/local/preprocess.nf'

workflow SPLINE {
    // Define a parameter to control preprocessing
    def skipPreprocessing = params.preprocessed

    // Conditionally run preprocessing
    if (!skipPreprocessing) {
        PREPROCESS(
            Channel.value(params.mmwrFile),
            Channel.value(params.projID)
        )
        clean_file = PREPROCESS.out.clean_csv
        use_preprocessed = true
    } else {
        clean_file = Channel.value(params.cleanFile ?: '')
        use_preprocessed = params.preprocessed
    }

    TRENDY(
      pathogens_ch,                          // 1. Pathogen channel (positional)
      Channel.value(params.mmwrFile),        // 2. mmwrFile (still needed for reference)
      Channel.value(params.censusFile_B),    // 3. censusFile_B
      Channel.value(params.censusFile_P),    // 4. censusFile_P
      Channel.value(params.travel),          // 5. travel
      Channel.value(params.cidt),            // 6. cidt
      Channel.value(params.projID),          // 7. projID
      Channel.value(params.whichScript),     // 8. whichScript (path to trendy.R)
      Channel.value(true),                   // 9. preprocessed flag (always true when we handle it here)
      clean_file                             // 10. cleanFile (path to cleaned CSV)
    )
}

