import tempfile
import os
import glob


# Config parameters from config.yaml file
TMP_DIR_ROOT = config["TMP_DIR_ROOT"]
RUN = config["RUN"]
SAMPLES_PE = config["samples_pe"] if "samples_pe" in config else []
SAMPLES_SE = config["samples_se"] if "samples_se" in config else []

# Path to programs
trimmomatic = config["software"]["trimmomatic"]
gzip        = config["software"]["gzip"]

# These rules are not processor intensive, and will execute on the head node
# without being allocated to compute nodes
localrules: raw_make_links_pe, raw_make_links_se, multiQC_run, multiQC_all

# This specifices the environment setup command for running compute jobs
# (for example, source activate kneaddata) and is specified in config.yaml
if "BOWTIE_ENV" in config:
    BOWTIE_ENV = config["BOWTIE_ENV"]
if "TRIM_ENV" in config:
    TRIM_ENV = config["TRIM_ENV"]
if "QC_ENV" in config:
    QC_ENV = config["QC_ENV"]

# DB info
if "HOST_DB" in config:
    HOST_DB = config["HOST_DB"]

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
            end = "R1 R2 U1 U2".split()
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

rule host_filter:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    """
    input:
        expand( # filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_PE,
                run = RUN,
                end = "R1 R2 U1 U2".split())


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
        set +u; {QC_ENV}; set -u
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
        unpaired_1 = "data/{sample}/{run}/trimmed/{sample}_U1.trimmed.fq.gz",
        unpaired_2 = "data/{sample}/{run}/trimmed/{sample}_U2.trimmed.fq.gz"
    params:
        forward  = "{sample}_R1.trimmed.fq.gz",
        reverse  = "{sample}_R2.trimmed.fq.gz",
        unpaired_1  = "{sample}_U1.trimmed.fq.gz",
        unpaired_2  = "{sample}_U2.trimmed.fq.gz",
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
                  set +u; {TRIM_ENV}; set -u
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
                  
                  scp %s/{params.forward} {output.forward}
                  scp %s/{params.reverse} {output.reverse}
                  scp %s/{params.unpaired_1} {output.unpaired_1}
                  scp %s/{params.unpaired_2} {output.unpaired_2}
                  """ % (temp_dir, temp_dir, temp_dir, temp_dir,
                         temp_dir, temp_dir, temp_dir, temp_dir
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
                  set +u; {TRIM_ENV}; set -u
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
        set +u; {QC_ENV}; set -u
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
        set +u; {QC_ENV}; set -u
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
         set +u; {QC_ENV}; set -u
         multiqc -f -o data/multiQC/all data 2> {log} 1>&2
         """


rule host_filter_pe:
    input:
        forward  = "data/{sample}/{run}/trimmed/{sample}_R1.trimmed.fq.gz",
        reverse  = "data/{sample}/{run}/trimmed/{sample}_R2.trimmed.fq.gz",
        unpaired_1 = "data/{sample}/{run}/trimmed/{sample}_U1.trimmed.fq.gz",
        unpaired_2 = "data/{sample}/{run}/trimmed/{sample}_U2.trimmed.fq.gz"
    output:
        forward  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        reverse  = "data/{sample}/{run}/host_filtered/{sample}_R2.trimmed.host_filtered.fq.gz",
        unpaired_1 = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        unpaired_2 = "data/{sample}/{run}/host_filtered/{sample}_U2.trimmed.host_filtered.fq.gz"
    params:
        forward_fn = "{sample}_R1.trimmed.host_filtered.fq",
        reverse_fn = "{sample}_R2.trimmed.host_filtered.fq",
        unpaired_1_fn = "{sample}_U1.trimmed.host_filtered.fq.gz",
        unpaired_2_fn = "{sample}_U2.trimmed.host_filtered.fq.gz"
    threads:
        8
    benchmark:
        "benchmarks/{run}/qc/host_filter_pe_{sample}.json"
    log:
        bowtie = "logs/{run}/bowtie/{sample}.log",
        other = "logs/{run}/qc/host_filter_pe_{sample}.log" 
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {BOWTIE_ENV}; set -u

                  bowtie2 -p {threads} -x {HOST_DB} --very-sensitive -1 {input.forward} -2 {input.reverse} 2> {log.bowtie}| \
                  samtools view -f 12 -F 256 2> {log.other}| \
                  samtools sort -@ {threads} -n 2> {log.other} | \
                  samtools view -bS 2> {log.other} | \
                  bedtools bamtofastq -i - -fq %s/{params.forward_fn} -fq2 %s/{params.reverse_fn} 2> {log.other}

                  {gzip} -c %s/{params.forward_fn} > {output.forward}
                  {gzip} -c %s/{params.reverse_fn} > {output.reverse}

                  bowtie2 -p {threads} -x {HOST_DB} --very-sensitive -U {input.unpaired_1} --un-gz %s/{params.unpaired_1_fn} -S /dev/null 2> {log.other}
                  bowtie2 -p {threads} -x {HOST_DB} --very-sensitive -U {input.unpaired_2} --un-gz %s/{params.unpaired_2_fn} -S /dev/null 2> {log.other}
                  scp %s/{params.unpaired_1_fn} {output.unpaired_1}
                  scp %s/{params.unpaired_2_fn} {output.unpaired_2}
                  """ % (temp_dir, temp_dir, temp_dir, temp_dir, temp_dir,
                         temp_dir, temp_dir, temp_dir))

