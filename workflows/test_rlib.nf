nextflow.enable.dsl = 2

process TestLibPaths {
    input:
      val dummy

    script:
    """
    echo "Printing R library paths and testing argparse..."
    Rscript -e 'print(".libPaths():"); print(.libPaths()); if(requireNamespace("argparse", quietly=TRUE)) { print("argparse loaded successfully") } else { print("argparse FAILED to load") }'
    """
}

workflow {
    TestLibPaths(Channel.value("test"))
}

