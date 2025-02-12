nextflow.enable.dsl = 2

// Create a channel with the list of pathogens to model.
// This has been modified to only run 2 pathogens for testing. It needs to have the others added before production.
Channel
    .from(['CAMPYLOBACTER', 'CYCLOSPORA'])
    .set { pathogens_ch }

// Include the TRENDY module from modules/local/trendy.nf
include { TRENDY } from '../modules/local/trendy.nf'

workflow SPLINE {
    TRENDY(
      pathogens_ch,                           // 1. Pathogen channel (positional)
      Channel.value(params.mmwrFile),         // 2. mmwrFile
      Channel.value(params.censusFile_B),     // 3. censusFile_B
      Channel.value(params.censusFile_P),     // 4. censusFile_P
      Channel.value(params.travel),           // 5. travel
      Channel.value(params.cidt),             // 6. cidt
      Channel.value(params.projID),           // 7. projID
      Channel.value(params.whichScript),      // 8. whichScript (path to trendy.R)
    // Channel.value(params.whichFunctions),   // 9. whichFunctions
      Channel.value(params.preprocessed),     // 10. preprocessed flag
      Channel.value(params.cleanFile ?: '')   // 11. cleanFile (path to cleaned CSV)
    )
}

