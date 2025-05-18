//
// This file holds several functions specific to the workflow/spline.nf in the CDCgov/FoodNetTrends pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowSpline {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        // Add parameter validation logic here if needed
    }

    //
    // Generate workflow summary for logging purposes
    //
    public static String paramsSummaryLog(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "  ${group}:\n"
                for (param in group_params.keySet()) {
                    summary_section += "    ${param}: ${group_params.get(param) ?: 'N/A'}\n"
                }
                summary_section += "\n"
            }
        }

        return summary_section
    }

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
    // Generate methods description for MultiQC
    //
    public static String toolCitationText(params) {
        // TODO nf-core: Optionally add tool citation text
        def citation_text = ""
        citation_text += "TODO Add tool citations here if available. Make sure to cite any databases or datasets that were used."
        return citation_text
    }

    public static String methodsDescriptionText(params) {
        // TODO nf-core: Optionally add methods text
        def methods_text = ""
        methods_text += "TODO Add methods description here if available."
        return methods_text
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "==================================================================================="
            System.exit(1)
        }
    }
}
