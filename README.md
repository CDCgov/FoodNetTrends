## Introduction

**cdc/spline** is a bioinformatics pipeline that performs spline modeling.

TODO: Write steps
1. Does thing1
2. Does thing2

## Usage
Run the pipeline by giving the following parameters:

    - outdir: absolute path to the output dir location
      - Example: /path/to/output/test240706
   - mmwrFile: absolute path to the MMWR file
      - Example: /path/to/mmwr9623_Jan2024.sas7bdat
	- censusFile_B: absolute path to the bacterial census file
      - Example: /path/to/cen9623.sas7bdat
	- censusFile_P: absolute path to the pathogenic census file
      - Example: /path/to/cen9623_para.sas7bdat
	- travel: list of travel status to include
      - Example: "NO,UNKNOWN,YES"
	- cidt: list of CIDT variables
      - Example: "CIDT+,CX+,PARASITIC"
	- projID: unique project identifier
      - Example: "20240706"

There are two ways to run the pipeline

1. Directly calls Nextflow
```bash
module load nextflow singularity conda

nextflow run main.nf \
	-entry CDC_SPLINE \
	-profile singularity,conda \
	-with-conda \
	-work-dir /path/to/output/work \
	--outdir /path/to/output \
	--mmwrFile "/path/to/mmwr9623_Jan2024.sas7bdat" \
	--censusFile_B "/path/to/cen9623.sas7bdat" \
	--censusFile_P "/path/to/cen9623_para.sas7bdat" \
	--travel "NO,UNKNOWN,YES" \
	--cidt "CIDT+,CX+,PARASITIC" \
	--projID "20240706"
```

2. Utilize the `run_workflow.sh` script.
   A. Update the `run_workflow.sh` script `outDir` and `dataDir` paths.
   B. Run the script
      ```bash
      bash run_workflow.sh run
      ```

## Credits

This pipeline was originally written by the SciComp .

`Spline` was largely developed by [OAMD's](https://www.cdc.gov/amd/index.html) SciComp Team, based on the scripts from [TBD](https://www.cdc.gov). Detailed contributions can be found in our [user-guides]().

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).