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
        mash = config['PARAMS']['MASH']
    threads:
        1
    log:
        "logs/{run}/analysis/mash_sketch_{sample}.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_sketch_{sample}.json"
    run:
        with tempfile.TemporaryDirectory(dir=TMP_DIR_ROOT) as temp_dir:
            shell("""
                  seqtk seq -a {input.forward} > {temp_dir}/R1.fna
                  seqtk seq -a {input.reverse} > {temp_dir}/R2.fna
                  seqtk seq -a {input.unpaired_1} > {temp_dir}/U1.fna
                  seqtk seq -a {input.unpaired_2} > {temp_dir}/U2.fna

                  cat {temp_dir}/R1.fna {temp_dir}/R2.fna {temp_dir}/U1.fna {temp_dir}/U2.fna > {temp_dir}/cat.fna

                  mash sketch {params.mash} {temp_dir}/cat.fna

                  mv {temp_dir}/cat.fna.msh {output}
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
        "combined_analysis/{run}/mash/mash.dist.txt"
    threads:
        1
    log:
        "logs/{run}/analysis/mash_dm.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_dm.json"
    run:
        for i in range(len(input)):
            for j in range(i,len(input)):
                shell("""
                      mash dist {input[i]} {input[j]} >> {output}
                      """)

rule mash_dm_write:
    """
    Writes square distance matrices from p values and distances that Mash makes
    """
    input:
        "combined_analysis/{run}/mash/mash.dist.txt"
    output:
        dist_matrix = "combined_analysis/{run}/mash/mash.dist.dm",
        p_matrix = "combined_analysis/{run}/mash/mash.dist.p"
    threads:
        1
    log:
        "logs/{run}/analysis/mash_dm_write.log"
    benchmark:
        "benchmarks/{run}/analysis/mash_dm_write.json"
    run:
        from skbio.stats.distance import DissimilarityMatrix
        import pandas as pd

        mash_vec = pd.read_csv(input[0], sep = '\t', header=None)

        # get sorted list of samples
        samples = sorted(set(mash_vec[0]) | set(mash_vec[1]))

        dm = np.zeros([len(samples),len(samples)])
        pm = np.zeros([len(samples),len(samples)])

        for s1, s2, d, p in zip(mash_vec[0],mash_vec[1],mash_vec[2],mash_vec[3]):
            i1 = samples.index(s1)
            i2 = samples.index(s2)
            print('s1: %s, s2: %s, i1: %s, i2: %s, d: %s, p: %s' % (s1, s2, i1, i2, d, p))
            dm[i1,i2] = d
            dm[i2,i1] = d
            pm[i1,i2] = p
            pm[i2,i1] = p

        sk_dm = DissimilarityMatrix(dm, ids=samples)
        sk_pm = DissimilarityMatrix(pm, ids=samples)

        sk_dm.write(output['dist_matrix'])
        sk_pm.write(output['p_matrix'])


