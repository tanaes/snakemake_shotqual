# how to install this workflow

# make a conda env for kneaddata
conda create -n kneaddata --yes python=3 pip pyyaml xlrd pandas jupyter notebook h5py

# install kneaddata and dependencies
source activate kneaddata
conda install -c bioconda --yes bowtie2 fastqc trimmomatic diamond samtools bedtools
pip install snakemake multiqc kneaddata biom-format

# download the kneaddata human genome database
mkdir -p ~/share/kd_dbs
kneaddata_database --download human bowtie2 ~/share/kd_dbs/human
cd ~/share

# Download metaphlan2
wget https://bitbucket.org/biobakery/metaphlan2/get/default.zip
unzip default.zip -d ~/share
mv ~/share/biobakery-metaphlan* ~/share/metaphlan2

# install humann2
pip install humann2

# download humann2 dbs
mkdir ~/share/humann2_db
humann2_databases --download chocophlan full ~/share/humann2_db
humann2_databases --download uniref uniref90_ec_filtered_diamond ~/share/humann2_db
