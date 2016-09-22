# Snakemake pipeline for shotgun data QC

A Snakemake workflow for assessing the quality of many shotgun sequence samples
simultaneously using [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)
and [MultiQC](http://multiqc.info).

It modifies some code from jlanga's snakemake khmer/trinity
[workflow](https://github.com/jlanga/khmer_trinity_snakemake). 


### Installation and usage

Currently, this tool depends on a recent version of Snakemake and runs:
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (v. 0.36)
- [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (≥ v. 0.11)
- [MultiQC](http://multiqc.info)

Snakemake reads the file config.yaml for several variables, including location
of relevant executables and parameters for Trimmomatic execution.

Forward and reverse reads currently have to be directly specified per sample
in the config.yaml file. For files produced by the UCSD IGM facility, I have
included a Jupyter Notebook to create this file from the sample manifest Excel
spreadsheet provided to IGM.

If executed directly using `snakemake --snakefile Snakefile --configfile config.yaml`,
the entire workflow will be run locally. Alternatively, the worklfow can be
executed in a Torque cluster environment with `bash launch.sh [dir] 
--snakefile Snakefile --configfile config.yaml`. This will parse cluster.json
for job submission parameters and run each rule instance as a separate cluster
job.


### Usage notes

Snakemake requires that inputs and outputs be specified explicitly to construct
the workflow graph. However, I have found that in practice there are samples
in the sequencing manifest that do not show up in the demultiplexed sequence
files, or do not yield all possible output files in Trimmomatic (for example,
if there are no R1 reads that survive trimming). For this reason, I have been
invoking this workflow with the `--keep-going` flag, which will run subsequent
steps even if not all outputs are successfully generated.