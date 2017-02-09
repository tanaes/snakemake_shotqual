import tempfile
import os
import glob


# Config parameters from config.yaml file
TMP_DIR_ROOT = config["TMP_DIR_ROOT"]
RUN = config["RUN"]
SAMPLES_PE = config["samples_pe"] if "samples_pe" in config else []
SAMPLES_SE = config["samples_se"] if "samples_se" in config else []

# Path to programs
trimmomatic = config["SOFTWARE"]["trimmomatic"]
gzip        = config["SOFTWARE"]["gzip"]

# These rules are not processor intensive, and will execute on the head node
# without being allocated to compute nodes
localrules: raw_make_links_pe, raw_make_links_se, multiQC_run, multiQC_all

# This specifices the environment setup command for running compute jobs
# (for example, source activate kneaddata) and is specified in config.yaml
if "BOWTIE_ENV" in config["ENVS"]:
    BOWTIE_ENV = config["ENVS"]["BOWTIE_ENV"]
if "TRIM_ENV" in config["ENVS"]:
    TRIM_ENV = config["ENVS"]["TRIM_ENV"]
if "QC_ENV" in config["ENVS"]:
    QC_ENV = config["ENVS"]["QC_ENV"]
if "HUMANN2_ENV" in config["ENVS"]:
    HUMANN2_ENV = config["ENVS"]["HUMANN2_ENV"]
if "METAPHLAN_ENV" in config["ENVS"]:
    METAPHLAN_ENV = config["ENVS"]["METAPHLAN_ENV"]
if "QUAST_ENV" in config["ENVS"]:
    QUAST_ENV = config["ENVS"]["QUAST_ENV"]

# Host DB specified in config file accessed in params section of rule

# Trimmomatic parameters accessed in params section of rule

# HUMAnN2 parameters accessed in params section of rule


#### Top-level rules: rules to execute a subset of the pipeline

rule all:
    """
    Rule to do all the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    And Host Filtering:
        - host_filter_pe
    And HUMANn2:
        - metaphlan2_sample_pe
        - combine_metaphlan
        - humann2_sample_pe
        - humann2_combine_tables
        - humann2_remove_unmapped
        - humann2_split_stratified_tables
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
        expand( # fastqc zipa nd html for trimmed PE data
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand( # fastqc zipa nd html for trimmed SE data
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        ),
        expand( # host-filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_PE,
                run = RUN,
                end = "R1 R2 U1 U2".split()) +
        expand( # filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_SE,
                run = RUN,
                end = "SE".split()),
        expand( # HUMANn2 on each sample
               "data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
               norm = config['PARAMS']['HUMANN2']['NORMS'],
               sample = SAMPLES_PE,
               run = RUN),
        expand( # Summarized HUMANn2 results
               "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_unstratified.biom",
               norm = config['PARAMS']['HUMANN2']['NORMS'],
               run = RUN,
               mapped=['all','mapped']),
        expand( # MultiQC for just this run
            "data/multiQC/{run}/multiqc_report.html",
            run = RUN
        ),
        "data/multiQC/all/multiqc_report.html"

rule basic_qc:
    """
    Rule to do basic qc, including host-filtering.

    All the Quality Control:
        - raw_fastqc
        - qc_trimmomatic_pe
        - qc_trimmomatic_se
        - qc_interleave_pe_pe
        - qc_fastqc
    And Host Filtering:
        - host_filter_pe
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
        expand( # fastqc zipa nd html for trimmed PE data
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_PE,
            end = "R1 R2".split(),
            run = RUN,
            extension = "zip html".split()
        ) + expand( # fastqc zipa nd html for trimmed SE data
            "data/{sample}/{run}/fastqc_trimmed/{sample}_{end}.trimmed_fastqc.{extension}",
            sample = SAMPLES_SE,
            end = "SE".split(),
            run = RUN,
            extension = "zip html".split()
        ),
        expand( # host-filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_PE,
                run = RUN,
                end = "R1 R2 U1 U2".split()) +
        expand( # filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_SE,
                run = RUN,
                end = "SE".split()),
        expand( # MultiQC for just this run
            "data/multiQC/{run}/multiqc_report.html",
            run = RUN
        ),
        "data/multiQC/all/multiqc_report.html"

rule host_filter:
    """
    Rule to do host-filtering
        - host_filter_pe
    """
    input:
        expand( # filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_PE,
                run = RUN,
                end = "R1 R2 U1 U2".split()) +
        expand( # filtered fastqs
                "data/{sample}/{run}/host_filtered/{sample}_{end}.trimmed.host_filtered.fq.gz",
                sample = SAMPLES_SE,
                run = RUN,
                end = "SE".split())

rule humann2:
    """
    Rule to do Humann2
        - metaphlan2_sample_pe
        - combine_metaphlan
        - humann2_sample_pe
        - humann2_combine_tables
        - humann2_remove_unmapped
        - humann2_split_stratified_tables
    """
    # params:
    #     norms = config['PARAMS']['HUMANN2']['NORMS']
    input:
        expand(# filtered fastqs
               "data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
               norm = config['PARAMS']['HUMANN2']['NORMS'],
               sample = SAMPLES_PE,
               run = RUN),
        expand(# stratified
               "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_unstratified.biom",
               norm = config['PARAMS']['HUMANN2']['NORMS'],
               run = RUN,
               mapped=['all','mapped'])



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
    """
    Makes symlinks from raw sequences files to the analysis directory.
    """
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
    """
    Makes symlinks from raw sequences files to the analysis directory.
    """
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
    """
    Do FASTQC reports for raw fastq sequence files. 
    One thread per fastq.gz file.
    """
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
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_raw {input.fastq} 2> {log} 1>&2
        """


rule qc_fastqc:
    """
    Do FASTQC reports for Trimmomatic-trimmed files. 
    One thread per fastq.gz file.
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
            --outdir data/{wildcards.sample}/{wildcards.run}/fastqc_trimmed {input.fastq} 2> {log} 1>&2
        """


rule qc_trimmomatic_pe:
    """
    Run trimmomatic on paired end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    
    Retains unpaired forward and reverse reads. 

    To avoid breaking, touches output file names in case there don't end up
    being any outputs in a particular category (sometimes happens with very
    shallowly sequenced samples). In this case, these filenames will exist but
    be empty.
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
        trimmomatic_params = config['PARAMS']['TRIMMOMATIC']["QUAL"],
        trimmomatic_clip = config['PARAMS']['TRIMMOMATIC']["ILLUMINACLIP"]
    benchmark:
        "benchmarks/{run}/qc/trimmomatic_pe_{sample}.json"
    log:
        "logs/{run}/qc/trimmomatic_pe_{sample}.log" 
    threads:
        6
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {TRIM_ENV}; set -u

                  {trimmomatic} PE \
                    -threads {threads} \
                    -{params.phred} \
                    {input.forward} \
                    {input.reverse} \
                    {temp_dir}/{params.forward} \
                    {temp_dir}/{params.unpaired_1} \
                    {temp_dir}/{params.reverse} \
                    {temp_dir}/{params.unpaired_2} \
                    ILLUMINACLIP:{params.adaptor}:{params.trimmomatic_clip} \
                    {params.trimmomatic_params} \
                  2> {log}

                  touch {temp_dir}/{params.forward}
                  touch {temp_dir}/{params.unpaired_1}
                  touch {temp_dir}/{params.reverse}
                  touch {temp_dir}/{params.unpaired_2}

                  scp {temp_dir}/{params.forward} {output.forward}
                  scp {temp_dir}/{params.reverse} {output.reverse}
                  scp {temp_dir}/{params.unpaired_1} {output.unpaired_1}
                  scp {temp_dir}/{params.unpaired_2} {output.unpaired_2}
                  """)



rule qc_trimmomatic_se:
    """
    Run trimmomatic on single end mode to eliminate Illumina adaptors and 
    remove low quality regions and reads.
    """
    input:
        single = "data/{sample}/{run}/raw/{sample}_SE.fq.gz"
    output:
        single = "data/{sample}/{run}/trimmed/{sample}_SE.trimmed.fq.gz"
    params:
        single = "{sample}_SE.trimmed.fq.gz",
        adaptor     = lambda wildcards: config["samples_se"][wildcards.sample]["adaptor"],
        phred       = lambda wildcards: config["samples_se"][wildcards.sample]["phred"],
        trimmomatic_params = config['PARAMS']['TRIMMOMATIC']["QUAL"],
        trimmomatic_clip = config['PARAMS']['TRIMMOMATIC']["ILLUMINACLIP"]
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
                      {temp_dir}/{params.single} \
                    ILLUMINACLIP:{params.adaptor}:{params.trimmomatic_clip} \
                    {params.trimmomatic_params} \
                  2> {log}

                  scp {temp_dir}/{params.single} {output.single} 
                  """)


rule multiQC_run:
    """
    Perform MultiQC summary on set of fastqc and logged output for the current
    analyzed run. 

    Will pass an explicit list of FastQC file outputs, in addition to the log
    directory. MultiQC will scan the log files and detect Trimmomatic and Bowtie
    outputs. 
    """
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
        multiqc -f -o data/multiQC/{wildcards.run} data/*/{wildcards.run} logs/{wildcards.run} 2> {log} 1>&2
        """


rule multiQC_all:
    """
    Perform MultiQC summary on set of fastqc and logged output for all runs
    in current project. 

    Passes the entire data directory, in addition to the log directory. MultiQC
    will scan the log files and detect Trimmomatic and Bowtie outputs, plus
    the FastQC outputs in the data directory. 
    """
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
         multiqc -f -o data/multiQC/all data logs 2> {log} 1>&2
         """


rule host_filter_pe:
    """
    Performs host read filtering on paired end data using Bowtie and Samtools/
    BEDtools. Takes the four output files generated by Trimmomatic. 

    Also requires an indexed reference (path specified in config). 

    First, uses Bowtie output piped through Samtools to only retain read pairs
    that are never mapped (either concordantly or just singly) to the indexed
    reference genome. Fastqs from this are gzipped into matched forward and 
    reverse pairs. 

    Unpaired forward and reverse reads are simply run through Bowtie and
    non-mapping gzipped reads output.

    All piped output first written to localscratch to avoid tying up filesystem.
    """
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
        unpaired_2_fn = "{sample}_U2.trimmed.host_filtered.fq.gz",
        host_db = config["HOST_DB"]
    threads:
        12
    benchmark:
        "benchmarks/{run}/qc/host_filter_pe_{sample}.json"
    log:
        bowtie = "logs/{run}/bowtie/{sample}.log",
        other = "logs/{run}/qc/host_filter_pe_{sample}.log" 
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {BOWTIE_ENV}; set -u

                  bowtie2 -p {threads} -x {params.host_db} --very-sensitive -1 {input.forward} -2 {input.reverse} 2> {log.bowtie}| \
                  samtools view -f 12 -F 256 2> {log.other}| \
                  samtools sort -T {temp_dir} -@ {threads} -n 2> {log.other} | \
                  samtools view -bS 2> {log.other} | \
                  bedtools bamtofastq -i - -fq {temp_dir}/{params.forward_fn} -fq2 {temp_dir}/{params.reverse_fn} 2> {log.other}

                  {gzip} -c {temp_dir}/{params.forward_fn} > {temp_dir}/{params.forward_fn}.gz
                  {gzip} -c {temp_dir}/{params.reverse_fn} > {temp_dir}/{params.reverse_fn}.gz

                  scp {temp_dir}/{params.forward_fn}.gz {output.forward}
                  scp {temp_dir}/{params.reverse_fn}.gz {output.reverse} 

                  bowtie2 -p {threads} -x {params.host_db} --very-sensitive -U {input.unpaired_1} --un-gz {temp_dir}/{params.unpaired_1_fn} -S /dev/null 2> {log.other}
                  bowtie2 -p {threads} -x {params.host_db} --very-sensitive -U {input.unpaired_2} --un-gz {temp_dir}/{params.unpaired_2_fn} -S /dev/null 2> {log.other}
                  scp {temp_dir}/{params.unpaired_1_fn} {output.unpaired_1}
                  scp {temp_dir}/{params.unpaired_2_fn} {output.unpaired_2}
                  """)



rule host_filter_se:
    """
    Performs host read filtering on single end data using Bowtie and Samtools/
    BEDtools. Takes the single output file generated by Trimmomatic. 

    Also requires an indexed reference (path specified in config). 

    Unpaired forward reads are simply run through Bowtie and
    non-mapping gzipped reads output.

    All piped output first written to localscratch to avoid tying up filesystem.
    """
    input:
        single = "data/{sample}/{run}/trimmed/{sample}_SE.trimmed.fq.gz"
    output:
        single  = "data/{sample}/{run}/host_filtered/{sample}_SE.trimmed.host_filtered.fq.gz"
    params:
        single = "{sample}_SE.trimmed.host_filtered.fq",
        host_db = config["HOST_DB"]
    threads:
        12
    benchmark:
        "benchmarks/{run}/qc/host_filter_pe_{sample}.json"
    log:
        bowtie = "logs/{run}/bowtie/{sample}.log",
        other = "logs/{run}/qc/host_filter_se_{sample}.log" 
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {BOWTIE_ENV}; set -u

                  bowtie2 -p {threads} -x {params.host_db} --very-sensitive -U {input.single} --un-gz {temp_dir}/{params.single} -S /dev/null 2> {log.other}
                  scp {temp_dir}/{params.single} {output.single}
                  """)


rule metaphlan2_sample_pe:
    """
    Runs MetaPhlan2 on a set of samples to create a joint taxonomic profile for
    input into HUMAnN2, based on the thinking that it is preferable to have a
    consistent Chocophlan reference database for the whole set of samples. This
    is especially true for shallowly sequenced samples. 

    Going to do just R1 reads for now. Because of how I've split PE vs SE
    processing and naming, still will need to make a separate rule for PE. 
    """
    input:
        paired_f  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        unpaired_f = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz"
    output:
        "data/{sample}/{run}/metaphlan2/{sample}_metaphlan_output.tsv"
    params:
        metaphlan_dir = config['PARAMS']['HUMANN2']["METAPHLAN_DIR"]
    threads:
        4
    log:
        "logs/{run}/analysis/metaphlan2_sample_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/metaphlan2_sample_pe_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {METAPHLAN_ENV}; set -u

                  zcat {input.paired_f} {input.unpaired_f} > {temp_dir}/input.fastq

                  {params.metaphlan_dir}/metaphlan2.py {temp_dir}/input.fastq \
                    --input_type fastq \
                    --mpa_pkl {params.metaphlan_dir}/db_v20/mpa_v20_m200.pkl \
                    --bowtie2db {params.metaphlan_dir}/db_v20/mpa_v20_m200 \
                    --nproc {threads} \
                    --tmp_dir {temp_dir} \
                    --no_map \
                    --input_type fastq > {output}  2> {log}
                  """)


rule combine_metaphlan:
    """
    Combines MetaPhlan2 output for unified taxonomic profile for Humann2.
    """
    input:
        expand("data/{sample}/{run}/metaphlan2/{sample}_metaphlan_output.tsv",
               sample=SAMPLES_PE, run=RUN)
    output:
        joint_prof = "data/combined_analysis/{run}/humann2/joined_taxonomic_profile.tsv",
        max_prof = "data/combined_analysis/{run}/humann2/joined_taxonomic_profile_max.tsv"
    threads:
        1
    log:
        "logs/{run}/analysis/combine_metaphlan.log"
    benchmark:
        "benchmarks/{run}/analysis/combine_metaphlan.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            for file in input:
                shell("cp {0} {1}/.".format(file, temp_dir))
            shell("""
                  set +u; {HUMANN2_ENV}; set -u

                  humann2_join_tables --input {temp_dir} --output {output.joint_prof} 2> {log} 1>&2
                  humann2_reduce_table --input {output.joint_prof} \
                  --output {output.max_prof} --function max --sort-by level 2>> {log} 1>&2
                  """)


rule humann2_sample_pe:
    """
    Runs HUMAnN2 pipeline using general defaults.

    Other HUMAnN2 parameters can be specified as a quoted string in 
    PARAMS: HUMANN2: OTHER. 

    Going to do just R1 reads for now. Because of how I've split PE vs SE
    processing and naming, still will need to make a separate rule for PE. 
    """
    input:
        paired_f  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        unpaired_f = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        metaphlan_in = "data/combined_analysis/{run}/humann2/joined_taxonomic_profile_max.tsv"
    output:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance.biom"
    params:
        metaphlan_dir = config['PARAMS']['HUMANN2']["METAPHLAN_DIR"],
        humann2_nt_db = config['PARAMS']['HUMANN2']["HUMANN2_NT_DB"],
        humann2_aa_db = config['PARAMS']['HUMANN2']["HUMANN2_AA_DB"],
        other = config['PARAMS']['HUMANN2']['OTHER']
    threads:
        8
    log:
        "logs/{run}/analysis/humann2_sample_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_sample_pe_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  set +u; {HUMANN2_ENV}; set -u

                  zcat {input.paired_f} {input.unpaired_f} > {temp_dir}/input.fastq

                  humann2 --input {temp_dir}/input.fastq \
                  --output {temp_dir}/{wildcards.sample} \
                  --output-basename {wildcards.sample} \
                  --nucleotide-database {params.humann2_nt_db} \
                  --protein-database {params.humann2_aa_db} \
                  --taxonomic-profile {input.metaphlan_in} \
                  --metaphlan {params.metaphlan_dir} \
                  --o-log {log} \
                  --threads {threads} \
                  --output-format biom {params.other} 2> {log} 1>&2


                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_genefamilies.biom {output.genefamilies}
                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_pathcoverage.biom {output.pathcoverage}
                  scp {temp_dir}/{wildcards.sample}/{wildcards.sample}_pathabundance.biom {output.pathabundance}
                  """)

rule humann2_renorm_tables:
    """
    Renormalizes HUMAnN2 per-sample tables, per recommendation in the HUMAnN2
    website. 

    Counts-per-million (cpm) or Relative Abundance (Relabund) can be specified
    as a list in the PARAMS: HUMANN2: NORMS variable in the config file.
    """
    input:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance.biom"
    output:
        genefamilies = "data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
        pathcoverage = "data/{sample}/{run}/humann2/{sample}_pathcoverage_{norm}.biom",
        pathabundance = "data/{sample}/{run}/humann2/{sample}_pathabundance_{norm}.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_renorm_tables_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_renorm_tables_{sample}.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u

              humann2_renorm_table --input {input.genefamilies} \
              --output {output.genefamilies} \
              --units {wildcards.norm} 2> {log} 1>&2

              humann2_renorm_table --input {input.pathcoverage} \
              --output {output.pathcoverage} \
              --units {wildcards.norm} 2>> {log} 1>&2

              humann2_renorm_table --input {input.pathabundance} \
              --output {output.pathabundance} \
              --units {wildcards.norm} 2>> {log} 1>&2
              """)


rule humann2_combine_tables:
    """
    Combines the per-sample normalized tables into a single run-wide table. 

    Because HUMAnN2 takes a directory as input, first copies all the individual
    tables generated in this run to a temp directory and runs on that.
    """
    input:
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_genefamilies_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm),
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_pathcoverage_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm),
        lambda wildcards: expand("data/{sample}/{run}/humann2/{sample}_pathabundance_{norm}.biom",
               sample=SAMPLES_PE, run=RUN, norm=wildcards.norm)
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_all.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_all.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_all.biom"
    log:
        "logs/{run}/analysis/humann2_combine_tables_{norm}.log"
    benchmark:
        "benchmarks/{run}/humann2_combine_tables_{norm}.log"
    run:
        with tempfile.TemporaryDirectory(dir='data/combined_analysis') as temp_dir:
            for file in input:
                shell("cp {0} {1}/.".format(file, temp_dir))
            shell("""
                  set +u; {HUMANN2_ENV}; set -u

                  humann2_join_tables --input {temp_dir} \
                  --output {output.genefamilies} \
                  --file_name genefamilies 2> {log} 1>&2

                  humann2_join_tables --input {temp_dir} \
                  --output {output.pathcoverage} \
                  --file_name pathcoverage 2>> {log} 1>&2

                  humann2_join_tables --input {temp_dir} \
                  --output {output.pathabundance} \
                  --file_name pathabundance 2>> {log} 1>&2

                  """)

rule humann2_remove_unmapped:
    """
    By default, HUMAnN2 includes the un-annoated reads (either unmapped in the
    first step of the pipeline or not matched in the translated alignment step)
    in the output files. In my experience, this causes relatively small
    differences in run quality (e.g. different read lengths) to have huge
    effects on the evaluated outcome, as the overall proportion of unmatched
    reads varies by run, and the compositionality of the data then causes large
    fluctuations in the count estimates of the annotated genes and pathways.

    To remove this problem, this rule renormalizes the combined tables after
    extracting unmatched read categories.
    """
    input:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_all.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_all.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_all.biom"
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_mapped.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_mapped.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_mapped.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_remove_un.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_remove_un.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u

              humann2_renorm_table --input {input.genefamilies} \
              --output {output.genefamilies} \
              --units {wildcards.norm} -s n 2> {log} 1>&2

              humann2_renorm_table --input {input.pathcoverage} \
              --output {output.pathcoverage} \
              --units {wildcards.norm} -s n 2>> {log} 1>&2

              humann2_renorm_table --input {input.pathabundance} \
              --output {output.pathabundance} \
              --units {wildcards.norm} -s n 2>> {log} 1>&2
              """)


rule humann2_split_stratified_tables:
    """
    Splits the grouped tables into separate grouped taxonomy-stratified and un-
    stratified tables for downstream analysis. Does this for both the combined
    tables including unmatched reads (*_all) and those excluding unmatched
    reads (*_mapped).

    The un-stratified tables should then be directly useful for downstream
    analysis in e.g. beta diversity. 
    """
    input:
        genefamilies = "data/combined_analysis/{run}/humann2/combined_genefamilies_{norm}_{mapped}.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/combined_pathcoverage_{norm}_{mapped}.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/combined_pathabundance_{norm}_{mapped}.biom"
    output:
        genefamilies = "data/combined_analysis/{run}/humann2/stratified/combined_genefamilies_{norm}_{mapped}_stratified.biom",
        pathcoverage = "data/combined_analysis/{run}/humann2/stratified/combined_pathcoverage_{norm}_{mapped}_stratified.biom",
        pathabundance = "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_stratified.biom",
        genefamilies_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_genefamilies_{norm}_{mapped}_unstratified.biom",
        pathcoverage_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_pathcoverage_{norm}_{mapped}_unstratified.biom",
        pathabundance_unstrat = "data/combined_analysis/{run}/humann2/stratified/combined_pathabundance_{norm}_{mapped}_unstratified.biom"
    threads:
        1
    log:
        "logs/{run}/analysis/humann2_split_stratified_tables.log"
    benchmark:
        "benchmarks/{run}/analysis/humann2_split_stratified_tables.json"
    run:
        shell("""
              set +u; {HUMANN2_ENV}; set -u

              humann2_split_stratified_table --input {input.genefamilies} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2> {log} 1>&2

              humann2_split_stratified_table --input {input.pathcoverage} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2>> {log} 1>&2

              humann2_split_stratified_table --input {input.pathabundance} \
              --output data/combined_analysis/{wildcards.run}/humann2/stratified 2>> {log} 1>&2
              """)


#### Mash rules
rule mash:
    input:
        expand('data/{sample}/{run}/mash/{sample}.fna.msh',
               sample=SAMPLES_PE, run=RUN),
        expand('data/{sample}/{run}/mash/{sample}.refseq.txt',
               sample=SAMPLES_PE, run=RUN),
        expand("data/combined_analysis/{run}/mash/mash.dist.dm", run=RUN),
        expand("data/combined_analysis/{run}/mash/mash.dist.p", run=RUN)



rule mash_sketch:
    """
    Sketches a trimmed and host-filtered fastq file. 
    
    There is almost no documentation for this tool, so it's problematic.

    Relevant parameters might be:
    -b : use bloom filtering on kmers to reduce impact of low-freq erros.
    -m N: use explicit depth filtering on kmers (bigger memory impact than bloom)
    """
    input:
        forward  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        reverse  = "data/{sample}/{run}/host_filtered/{sample}_R2.trimmed.host_filtered.fq.gz",
        unpaired_1 = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        unpaired_2 = "data/{sample}/{run}/host_filtered/{sample}_U2.trimmed.host_filtered.fq.gz"
    output:
        'data/{sample}/{run}/mash/{sample}.fna.msh'
    params:
        mash = config['SOFTWARE']['mash'],
        seqtk = config['SOFTWARE']['seqtk'],
        mash_params = config['PARAMS']['MASH']['OTHER'],
        output_base = 'data/{sample}/{run}/mash/{sample}.fna'
    threads:
        1
    log:
        "logs/{run}/analysis/mash_sketch_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_sketch_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  {params.seqtk} seq -a {input.forward} > {temp_dir}/R1.fna
                  {params.seqtk} seq -a {input.reverse} > {temp_dir}/R2.fna
                  {params.seqtk} seq -a {input.unpaired_1} > {temp_dir}/U1.fna
                  {params.seqtk} seq -a {input.unpaired_2} > {temp_dir}/U2.fna

                  cat {temp_dir}/R1.fna {temp_dir}/R2.fna {temp_dir}/U1.fna {temp_dir}/U2.fna > {temp_dir}/{wildcards.sample}

                  {params.mash} sketch {params.mash_params} -o {params.output_base} {temp_dir}/{wildcards.sample}
                  """)


rule mash_refseq:
    """
    Compares a mash sketch against refseq sketch. 

    Requires that the sketches have same -k values -- for RefSeqDefault, 
    -k should equal 21. 
    """
    input:
        'data/{sample}/{run}/mash/{sample}.fna.msh'
    output:
        'data/{sample}/{run}/mash/{sample}.refseq.txt'
    params:
        mash = config['SOFTWARE']['mash'],
        db = config['PARAMS']['MASH']['REFSEQ_DB']
    threads:
        1
    log:
        "logs/{run}/analysis/mash_refseq_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_refseq_{sample}.json"
    run:
        shell("""
              {params.mash} dist {params.db} {input} | sort -gk3 > {output}
              """)

rule mash_dm:
    """
    Makes mash distance output file
    """
    input:
        expand( # trimmomatic output for PE data
            "data/{sample}/{run}/mash/{sample}.fna.msh",
            sample = SAMPLES_PE,
            run = RUN)
    output:
        "data/combined_analysis/{run}/mash/mash.dist.txt"
    params:
        mash = config['SOFTWARE']['mash']
    threads:
        1
    log:
        "logs/{run}/analysis/mash_dm.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_dm.json"
    run:
        for i in range(len(input)):
            for j in range(i,len(input)):
                thing1 = input[i]
                thing2 = input[j]
                shell("""
                      {params.mash} dist {thing1} {thing2} >> {output}
                      """)

rule mash_dm_write:
    """
    Writes square distance matrices from p values and distances that Mash makes
    """
    input:
        "data/combined_analysis/{run}/mash/mash.dist.txt"
    output:
        dist_matrix = "data/combined_analysis/{run}/mash/mash.dist.dm",
        p_matrix = "data/combined_analysis/{run}/mash/mash.dist.p"
    threads:
        1
    log:
        "logs/{run}/analysis/mash_dm_write.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_dm_write.json"
    run:
        from skbio.stats.distance import DissimilarityMatrix
        import pandas as pd
        import numpy as np

        mash_vec = pd.read_csv(input[0], sep = '\t', header=None)

        # get sorted list of samples
        samples = sorted(set(mash_vec[0]) | set(mash_vec[1]))

        dm = np.zeros([len(samples),len(samples)])
        pm = np.zeros([len(samples),len(samples)])

        # fill matrices with values
        for s1, s2, d, p in zip(mash_vec[0],mash_vec[1],mash_vec[2],mash_vec[3]):
            i1 = samples.index(s1)
            i2 = samples.index(s2)
            print('s1: %s, s2: %s, i1: %s, i2: %s, d: %s, p: %s' % (s1, s2, i1, i2, d, p))
            dm[i1,i2] = d
            dm[i2,i1] = d
            pm[i1,i2] = p
            pm[i2,i1] = p

        ids = [os.path.basename(x) for x in samples]
        sk_dm = DissimilarityMatrix(dm, ids=ids)
        sk_pm = DissimilarityMatrix(pm, ids=ids)

        sk_dm.write(output['dist_matrix'])
        sk_pm.write(output['p_matrix'])


rule assemble:
    input:
        expand('data/{sample}/{run}/megahit/{sample}.contigs.fa',
               sample=SAMPLES_PE, run=RUN),
        expand('data/{sample}/{run}/metaquast/done.txt',
               sample=SAMPLES_PE, run=RUN),
        expand("data/{sample}/{run}/quast/done.txt",
               sample=SAMPLES_PE, run=RUN)



rule megahit:
    """
    Run Megahit assembly on fastq
    """
    input:
        forward  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        reverse  = "data/{sample}/{run}/host_filtered/{sample}_R2.trimmed.host_filtered.fq.gz",
        unpaired_1 = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        unpaired_2 = "data/{sample}/{run}/host_filtered/{sample}_U2.trimmed.host_filtered.fq.gz"
    output:
        'data/{sample}/{run}/megahit/{sample}.contigs.fa'
    params:
        megahit = config['SOFTWARE']['megahit'],
        memory = 128
    threads:
        8
    log:
        "logs/{run}/assembly/megahit_{sample}.log"
    benchmark:
        "benchmarks/{run}/assembly/megahit_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            mem_b = params.memory * 1000000000
            outdir = os.path.dirname(output[0])
            shell("""
                  rm -r {outdir}

                  {params.megahit} -1 {input.forward} -2 {input.reverse} -r {input.unpaired_1} -r {input.unpaired_2} \
                  -m {mem_b} -t {threads} -o {outdir} --out-prefix {wildcards.sample} --tmp-dir {temp_dir} 2> {log} 1>&2
                  """)


rule metaquast:
    input:
        'data/{sample}/{run}/megahit/{sample}.contigs.fa'
    output:
        "data/{sample}/{run}/metaquast/done.txt"
    log:
        "logs/{run}/assembly/metaquast_{sample}.log"
    params:
        metaquast = config['SOFTWARE']['metaquast']
    threads:
        8
    run:
        outdir = os.path.dirname(output[0])
        shell("""
              set +u; {QUAST_ENV}; set -u

              {params.metaquast} -t {threads} -o {outdir} {input} 2> {log} 1>&2

              touch {output}
              """)


rule quast:
    input:
        'data/{sample}/{run}/megahit/{sample}.contigs.fa'
    output:
        "data/{sample}/{run}/quast/done.txt"
    log:
        "logs/{run}/assembly/quast_{sample}.log"
    params:
        quast = config['SOFTWARE']['quast']
    threads:
        8
    run:
        outdir = os.path.dirname(output[0])
        shell("""
              set +u; {QUAST_ENV}; set -u

              {params.quast} -t {threads} -o {outdir} {input} 2> {log} 1>&2

              touch {output}
              """)






