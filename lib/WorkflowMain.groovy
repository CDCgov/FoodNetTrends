//
// This file holds several functions specific to the main.nf workflow in the CDCgov/FoodNetTrends pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://github.com/CDCgov/FoodNetTrends\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/CDCgov/FoodNetTrends/blob/master/CITATIONS.md"
    }

    //
    // Print help to screen if required
    //
    public static void help(workflow, params) {
        // Don't need to do anything here since we handle help in main.nf directly
    }

    //
    // Print parameter summary log to screen
    //
    public static void summary(workflow, params, log) {
        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = workflow.manifest.version ?: "DEV"
            log.info "${workflow.manifest.name} v${workflow_version}"
            System.exit(0)
        }
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        help(workflow, params)

        // Print workflow version and exit on version parameter
        if (params.version) {
            String workflow_version = workflow.manifest.version ?: "DEV"
            log.info "${workflow.manifest.name} v${workflow_version}"
            System.exit(0)
        }

        // Print parameter summary
        summary(workflow, params, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        Utils.awsBatch(workflow, params)

        // Ensure workflow supports singularity if it's enabled
        Utils.checkSingularity(workflow, params)
    }
}
