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
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

The current version also has some provisions for executing [HUMAnN2](https://bitbucket.org/biobakery/humann2/wiki/Home)
for looking at functional profiles from the metagenomes. 

Snakemake reads the file config.yaml for several variables, including location
of relevant executables and parameters for Trimmomatic execution.

Forward and reverse reads currently have to be directly specified per sample
in the config.yaml file. For files produced by the UCSD IGM facility, I have
included a Jupyter Notebook to create this file from the sample manifest Excel
spreadsheet provided to IGM.

If executed directly using `snakemake --snakefile Snakefile --configfile config.yaml`,
the entire workflow will be run locally. Alternatively, the workflow can be
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

Currently, the test data are completely artificial very short sequences. So that
these run, the template config.yaml file has a very permissive length filter for
Trimmomatic (`MINLEN:3`). For real data, a better default is `MINLEN:32`.

Finally, note that disk access-intensive steps are set to run on a temporary
directory to allow execution on local scratch space in a cluster environment.
This variable is called `TMP_DIR_ROOT` in the config.yaml, and should be set to
the local scratch directory to enable this behavior. 


### How to run

On Barnacle, we want to avoid running compute-intensive jobs on the login node.
That's what happens if we just run the included Snakefile without any
additional information about how to access the cluster.

local execution (**DON'T DO THIS**):
```bash
snakemake --configfile config_Run1.yaml
```

Instead, I've provided a launch.sh script that is set up with some defaults
chosen to improve execution on our cluster. Here's how you run it:

cluster execution (**DO THIS**):
```bash
bash launch.sh ./ --configfile config_Run1.yaml
```

Here's what's goingon behind the scenes in launch.sh to invoke the Snakemake
workflow:

```bash
snakemake -j 16 \
--local-cores 4 \
-w 90 \
--max-jobs-per-second 8 \
--cluster-config cluster.json \
--cluster "qsub -k eo -m n -l nodes=1:ppn={cluster.n} -l mem={cluster.mem}gb -l walltime={cluster.time}" \
--directory "$@"
```

Let's go through what each of these parameters does.

`-j 16`: Runs no more than 16 jobs concurrently. If you have 96 samples that
each need to get FastQC'd, it will only run 16 of these jobs at a time.

`--local-cores 4`: For rules specified as local rules (like linking files),
limits to use of 4 CPUs at a time. 

`-w 90`: Waits for at most 90 seconds after a job executes for the output files
to be available. This has to do with tolerating latency on the filesystem:
sometimes a file is created by a job but isn't immediately visible to the
Snakemake process that's scheduling things.

`--max-jobs-per-second 8`: Limits the rate at which Snakemake is sending jobs
to the cluster.

`--cluster-config cluster.json` Looks in the current directory for a file
called cluster.json that contains information about how many resources to
request from the cluster for each rule type.

`--cluster "qsub -k eo [...]"`: This tells Snakemake how to send a job to the
cluster scheduler, and how to request the specific resources defined in the
cluster.json file.

`--directory "$@"`: This passes all the input provided after `bash launch.sh`
as further input to Snakemake. Because it comes right after the `--directory`
flag, it's going to expect the first element of that input to be the path to
the working directory where Snakemake should execute. 

