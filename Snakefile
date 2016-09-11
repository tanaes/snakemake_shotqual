import tempfile
import os
import glob

#configfile: "config.yaml.template"


RUNS = config["RUNS"]
SAMPLES = config["SAMPLES"]
TRIMMOMATIC_JAR = config["TRIMMOMATIC_JAR"]
FASTQC_PATH = config["FASTQC_PATH"]
TMP_DIR_ROOT = config["TMP_DIR_ROOT"]


def get_r1(wildcards):
    return glob.glob(os.path.join(RUNS[wildcards.run], SAMPLES[wildcards.sample] + "_*_R1_*.fastq.gz"))

def get_r2(wildcards):
    return glob.glob(os.path.join(RUNS[wildcards.run], SAMPLES[wildcards.sample] + "_*_R2_*.fastq.gz"))


#localrules: all, gather_assemblies

rule all:
    input:
        #expand("data/{sample}/{run}/raw/{sample}_R1.fq.gz", sample=SAMPLES.keys(), run=RUNS.keys()),
        #expand("data/{sample}/{run}/raw/{sample}_R2.fq.gz", sample=SAMPLES.keys(), run=RUNS.keys()),
        expand("data/{sample}/{run}/raw/FastQC/{sample}_R2_fastqc.html",  sample=SAMPLES.keys(), run=RUNS.keys()),
        expand("data/{sample}/{run}/trimmed/FastQC/{sample}_R1_paired_fastqc.html",  sample=SAMPLES.keys(), run=RUNS.keys()),
        #expand("data/multiQC/{run}/multiqc_report.html", run=RUNS.keys())

rule link_files:
    input: 
        r1 = get_r1,
        r2 = get_r2
    output: 
        r1 = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        r2 = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    shell:
        "ln -s {input.r1} {output.r1}; ln -s {input.r2} {output.r2}"

rule raw_fastqc:
    input:
        r1 = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        r2 = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    output:
        "data/{sample}/{run}/raw/FastQC/{sample}_R1_fastqc.zip",
        "data/{sample}/{run}/raw/FastQC/{sample}_R2_fastqc.zip",
        "data/{sample}/{run}/raw/FastQC/{sample}_R1_fastqc.html",
        "data/{sample}/{run}/raw/FastQC/{sample}_R2_fastqc.html"
    params:
        fastqc_path = FASTQC_PATH
    run:
        shell("""
              {params.fastqc_path} --outdir data/{wildcards.sample}/{wildcards.run}/raw/FastQC/ {input.r1} {input.r2}
              """)


## use trimmomatic to trim low quality bases and adaptors
rule clean_fastq:
    input:
        r1 = "data/{sample}/{run}/raw/{sample}_R1.fq.gz",
        r2 = "data/{sample}/{run}/raw/{sample}_R2.fq.gz"
    output:
        r1_p = "data/{sample}/{run}/trimmed/{sample}_R1_paired.fq.gz",
        r2_p = "data/{sample}/{run}/trimmed/{sample}_R2_paired.fq.gz"
    params:
        r1_u = "data/{sample}/{run}/trimmed/{sample}_R1_unpaired.fq.gz",
        r2_u = "data/{sample}/{run}/trimmed/{sample}_R2_unpaired.fq.gz",
        adapter = config["ADAPTER"],
        leading = config["LEADING"],
        trailing = config["TRAILING"],
        window = config["WINDOW"],
        minlen = config["MINLEN"],
        threads = 8,
        trimmomatic_jar = TRIMMOMATIC_JAR
    run:
        shell("""
              java -jar {params.trimmomatic_jar} PE -threads {params.threads} \
              {input.r1} {input.r2} \
              {output.r1_p} {params.r1_u} {output.r2_p} {params.r1_u} \
              ILLUMINACLIP:{params.adapter} LEADING:{params.leading} \
              TRAILING:{params.trailing} SLIDINGWINDOW:{params.window} \
              MINLEN:{params.minlen}
              """)


rule trimmed_fastqc:
    input:
        r1 = "data/{sample}/{run}/trimmed/{sample}_R1_paired.fq.gz",
        r2 = "data/{sample}/{run}/trimmed/{sample}_R2_paired.fq.gz"
    output:
        "data/{sample}/{run}/trimmed/FastQC/{sample}_R1_paired_fastqc.zip",
        "data/{sample}/{run}/trimmed/FastQC/{sample}_R2_paired_fastqc.zip",
        "data/{sample}/{run}/trimmed/FastQC/{sample}_R1_paired_fastqc.html",
        "data/{sample}/{run}/trimmed/FastQC/{sample}_R2_paired_fastqc.html"
    params:
        fastqc_path = FASTQC_PATH
    run:
        shell("""
              {params.fastqc_path} --outdir data/{wildcards.sample}/{wildcards.run}/trimmed/FastQC/ {input.r1} {input.r2}
              """)

rule multiQC_run:
    input: 
        lambda wildcards: expand("data/{sample}/{run}/trimmed/FastQC/{sample}_R1_paired_fastqc.zip", sample=SAMPLES.keys(), run=wildcards.run),
        lambda wildcards: expand("data/{sample}/{run}/trimmed/FastQC/{sample}_R2_paired_fastqc.zip", sample=SAMPLES.keys(), run=wildcards.run),
        lambda wildcards: expand("data/{sample}/{run}/raw/FastQC/{sample}_R1_fastqc.zip", sample=SAMPLES.keys(), run=wildcards.run),
        lambda wildcards: expand("data/{sample}/{run}/raw/FastQC/{sample}_R2_fastqc.zip", sample=SAMPLES.keys(), run=wildcards.run)
    output:
        "data/multiQC/{run}/multiqc_report.html"
    params:
        error = "data/errors.txt"
    shell:
        "source activate multiqc;multiqc -f -o data/multiQC/{wildcards.run} data/*/{wildcards.run} 2> errors.txt"

rule multiQC_all:
    input:
        expand("data/{sample}/{run}/trimmed/FastQC/{sample}_R1_paired_fastqc.zip", sample=SAMPLES.keys(), run=RUNS.keys()),
        expand("data/{sample}/{run}/trimmed/FastQC/{sample}_R2_paired_fastqc.zip", sample=SAMPLES.keys(), run=RUNS.keys()),
        expand("data/{sample}/{run}/raw/FastQC/{sample}_R1_fastqc.zip", sample=SAMPLES.keys(), run=RUNS.keys()),
        expand("data/{sample}/{run}/raw/FastQC/{sample}_R2_fastqc.zip", sample=SAMPLES.keys(), run=RUNS.keys())
    output:
        "data/multiQC/all/multiqc_report.html"
    shell:
        """
        source activate multiqc \
        multiqc data -o data/multiQC/all
        """
