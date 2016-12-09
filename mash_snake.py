
rule mash_sketch:
    input:
        forward  = "data/{sample}/{run}/host_filtered/{sample}_R1.trimmed.host_filtered.fq.gz",
        reverse  = "data/{sample}/{run}/host_filtered/{sample}_R2.trimmed.host_filtered.fq.gz",
        unpaired_1 = "data/{sample}/{run}/host_filtered/{sample}_U1.trimmed.host_filtered.fq.gz",
        unpaired_2 = "data/{sample}/{run}/host_filtered/{sample}_U2.trimmed.host_filtered.fq.gz"
    output:
        'data/{sample}/{run}/mash/{sample}.fna.msh'
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  seqtk seq -a {input.forward} > {temp_dir}/R1.fna
                  seqtk seq -a {input.reverse} > {temp_dir}/R2.fna
                  seqtk seq -a {input.unpaired_1} > {temp_dir}/U1.fna
                  seqtk seq -a {input.unpaired_2} > {temp_dir}/U2.fna

                  cat {temp_dir}/R1.fna {temp_dir}/R2.fna {temp_dir}/U1.fna {temp_dir}/U2.fna > {temp_dir}/cat.fna

                  mash sketch {temp_dir}/cat.fna

                  mv {temp_dir}/cat.fna.msh {output}
                  """)

rule mash_dm:
    """
    Makes mash distance matrix
    """
    input:
        expand( # trimmomatic output for PE data
            "data/{sample}/{run}/mash/{sample}.fna.msh",
            sample = SAMPLES_PE,
            run = RUN)
    output:
        "combined_analysis/{run}/mash/mash.dist.txt"
    threads:
        1
    log:
        "logs/{run}/raw/make_links_pe_{sample}.log"
    benchmark:
        "benchmarks/{run}/raw/make_links_pe_{sample}.json"
    run:
        for i in range(len(input)):
            for j in range(i,len(input)):
                shell("""
                      mash dist {input[i]} {input[j]} >> {output}
                      """)






