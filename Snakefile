import tempfile
import os
import glob


configfile: "config.yaml"

TMP_DIR_ROOT = config["TMP_DIR_ROOT"]
RUN = config["RUN"]
SAMPLES_PE = config["samples_pe"] if config["samples_pe"] is not None else []
SAMPLES_SE = config["samples_se"] if config["samples_se"] is not None else []

# Path to programs (or element on path)
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]


localrules: raw_make_links_pe, raw_make_links_se, multiQC_run, multiQC_all


#### Top-level rules: rules to execute a subset of the pipeline

rule all:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    """
    input:
        expand( # fastqc zip and html for raw PE data
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            run = RUN,
            end = "R1 R2".split(),
            extension = "zip html".split()
        ) + expand( # fastqc zip and html for raw SE data
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_SE,
            run = RUN,
            end = "SE".split(),
            extension = "zip html".split()
        ),
        expand( # trimmomatic output for PE data
            "data/{sample}/{run}/trimmed/{sample}_{end}.trimmed.fq.gz",
            sample = SAMPLES_PE,
            run = RUN,
            end = "R1 R2 up".split()
        ) + expand( # fastqc zip and html for raw SE data
            "data/{sample}/{run}/trimmed/{sample}_{end}.trimmed.fq.gz",
            sample = SAMPLES_SE,
            run = RUN,
            end = "SE".split()
        ),
        expand(
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand(
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        ),
        expand(
            "data/multiQC/{run}/multiqc_report.html",
            run = RUN
        ),
        "data/multiQC/all/multiqc_report.html"


rule raw_fastqc:
    """
    Rule to do just QC on raw reads:
        - raw_fastqc
    """
    input:
        expand(
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand(
            "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        )


rule raw_make_links_pe:
    input:
        forward = lambda wildcards: config["samples_pe"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples_pe"][wildcards.sample]["reverse"]
    output:
        forward = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        reverse = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    threads:
        1
    log:
        "logs/{run}/raw/make_links_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/raw/make_links_pe_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.forward}) {output.forward} 2> {log}
        ln -s $(readlink -f {input.reverse}) {output.reverse} 2>> {log}
        """ 


rule raw_make_links_se:
    input:
        single = lambda wildcards: config["samples_se"][wildcards.sample]["forward"],
    output:
        single = "data/{sample}/{run}/raw/{sample}_SE.fq.gz"
    threads:
        1
    log:
        "logs/{run}/raw/make_links_se_{sample}.log"
    benchmark:
        "benchmarks/{run}/raw/make_links_se_{sample}.json"
    shell:
        """
        ln -s $(readlink -f {input.single}) {output.single} 2>  {log}
        """ 


rule raw_fastqc_sample:
    input:
        fastq = "data/{sample}/{run}/raw/{sample}_{end}.fq.gz"
    output:
        html = "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip"
    threads:
        1
    params:
        html = "data/{sample}/{run}/raw/{sample}_{end}_fastqc.html",
        zip =  "data/{sample}/{run}/raw/{sample}_{end}_fastqc.zip"
    log:
        "logs/{run}/raw/fastqc_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/raw/fastqc_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_raw \
            {input.fastq} \
        2> {log} 1>&2
        """


rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    Inputs _1 and _2 are piped through gzip/pigz.
    Outputs _1 and _2 are piped to gzip/pigz (level 9).
    Outputs _3 and _4 are compressed with the builtin compressor from 
    Trimmomatic. Further on they are catted and compressed with gzip/pigz 
    (level 9).
    Note: The cut -f 1 -d " " is to remove additional fields in the FASTQ
    header. It is done posterior to the trimming since the output comes 
    slower than the input is read.
    Number of threads used:
        4 for trimmomatic
        2 for gzip inputs
        2 for gzip outputs
        Total: 8
    """
    input:
        forward = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        reverse = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    output:
        forward  = "data/{sample}/{run}/trimmed/{sample}_R1.trimmed.fq.gz",
        reverse  = "data/{sample}/{run}/trimmed/{sample}_R2.trimmed.fq.gz",
        unpaired = "data/{sample}/{run}/trimmed/{sample}_up.trimmed.fq.gz"
    params:
        forward  = "{sample}_R1.trimmed.fq.gz",
        reverse  = "{sample}_R2.trimmed.fq.gz",
        unpaired_1  = "{sample}_up_R1.trimmed.fq.gz",
        unpaired_2  = "{sample}_up_R2.trimmed.fq.gz",
        unpaired = "{sample}_up.trimmed.fq.gz",
        adaptor     = lambda wildcards: config["samples_pe"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_pe"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/{run}/qc/trimmomatic_pe_{sample}.json"
    log:
        "logs/{run}/qc/trimmomatic_pe_{sample}.log" 
    threads:
        8
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  {trimmomatic} PE \
                    -threads {threads} \
                    -{params.phred} \
                    {input.forward} \
                    {input.reverse} \
                    %s/{params.forward} \
                    %s/{params.unpaired_1} \
                    %s/{params.reverse} \
                    %s/{params.unpaired_2} \
                    ILLUMINACLIP:{params.adaptor}:2:30:10 \
                    {params.trimmomatic_params} \
                  2> {log}
                  
                  zcat %s/{params.unpaired_1} %s/{params.unpaired_2} |
                  cut -f 1 -d " " |
                  {gzip} -9 > %s/{params.unpaired}
                  
                  scp %s/{params.forward} {output.forward}
                  scp %s/{params.reverse} {output.reverse}
                  scp %s/{params.unpaired_2} {output.unpaired}
                  """ % (temp_dir, temp_dir, temp_dir, temp_dir, temp_dir,
                         temp_dir, temp_dir, temp_dir, temp_dir, temp_dir
                         ))



rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and 
        remove low quality regions and reads.
    Input is piped through gzip/pigz.
    Output is piped to gzip.
    Threads used:
        4 for trimmomatic
        1 for gzip input
        1 for gzip output
    """
    input:
        single = "data/{sample}/{run}/raw/{sample}_SE.fq.gz"
    output:
        single = "data/{sample}/{run}/trimmed/{sample}_SE.trimmed.fq.gz"
    params:
        single = "{sample}_SE.trimmed.fq.gz",
        adaptor     = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config["trimmomatic_params"]
    benchmark:
        "benchmarks/{run}/qc/trimmomatic_se_{sample}.json"
    log:
        "logs/{run}/qc/trimmomatic_se_{sample}.log" 
    threads:
        6
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  {trimmomatic} SE \
                      -threads {threads} \
                      -{params.phred} \
                      {input.single} \
                      %s/{params.single} \
                      ILLUMINACLIP:{params.adaptor}:2:30:10 \
                      {params.trimmomatic_params} \
                  2> {log}

                  scp %s/{params.single} {output.single} 
                  """ % (temp_dir, temp_dir))


rule qc_fastqc:
    """
    Do FASTQC reports
    One thread per fastq.gz file
    """
    input:
        fastq = "data/{sample}/{run}/trimmed/{sample}_{end}.trimmed.fq.gz"
    output:
        html = "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.html",
        zip =  "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.zip"
    threads:
        1
    log:
        "logs/{run}/qc/fastqc_trimmed_{sample}_{end}.log"
    benchmark:
        "benchmarks/{run}/qc/fastqc_trimmed_{sample}_{end}.json"
    shell:
        """
        fastqc \
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_trimmed \
            {input.fastq} 
        2> {log} 1>&2
        """


rule multiQC_run:
    input: 
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split())
    output:
        "data/multiQC/{run}/multiqc_report.html"
    threads:
        1
    log:
        "logs/{run}/qc/multiqc.log"
    benchmark:
        "benchmarks/{run}/qc/multiqc.json"
    shell:
        """
        set +u; source activate multiqc; set -u
        multiqc -f -o data/multiQC/{wildcards.run} data/*/{wildcards.run} 2> {log} 1>&2
        """


rule multiQC_all:
    input:
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_PE, run=RUN, end="R1 R2".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.html", sample=SAMPLES_SE, run=RUN, end="SE".split()),
        expand("data/{sample}/{run}/fastqc_raw/{sample}_{end}_fastqc.zip", sample=SAMPLES_SE, run=RUN, end="SE".split())
    output:
        "data/multiQC/all/multiqc_report.html"
    threads:
        1
    log:
        "logs/multiqc_all.log"
    benchmark:
        "benchmarks/multiqc_all.json"
    shell:
         """
         set +u; source activate multiqc; set -u
         multiqc -f -o data/multiQC/all data 2> {log} 1>&2
         """
