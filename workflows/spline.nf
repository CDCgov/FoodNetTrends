// Include the TRENDY module
include { TRENDY } from '../modules/local/trendy.nf'

workflow SPLINE {
    TRENDY(
        params.whichScript,
        params.whichFunctions,
        params.mmwrFile,
        params.censusFile_B,
        params.censusFile_P,
        params.travel,
        params.cidt,
        params.projID
    )
}

