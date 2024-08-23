# Mini-Project: RNA-seq Analysis - From Raw Data to Count Matrix 
# Create the main working directory
mkdir RNAseq_Project

# Navigate into the working directory
cd RNAseq_Project

# Create subdirectories for raw data, reference files, and results
mkdir RawData Reference Results 

# Download Data files and unzip them ( skipped )
wget http://example.com/path/to/

# install Unzip package tool
sudo apt-get install unzip
unzip Sequencing data.zip

# Create 2 files and Move FASTQ files into RawData directory
cd RNAseq/RAWDATA cd Sequencing data
mkdir Cancer Normal
mv cancer_sample_1.read1.fastq  RAWDATA/Sequencing data/cancer
mv cancer_sample_1.read2.fastq  RAWDATA/Sequencing data/cancer
mv cancer_sample_2.read1.fastq  RAWDATA/Sequencing data/cancer
mv cancer_sample_2.read2.fastq  RAWDATA/Sequencing data/cancer
mv cancer_sample_3.read1.fastq  RAWDATA/Sequencing data/cancer
mv cancer_sample_3.read2.fastq  RAWDATA/Sequencing data/cancer
mv normal_sample_1.read1.fastq  RAWDATA/Sequencing data/Normal
mv normal_sample_1.read2.fastq  RAWDATA/Sequencing data/Normal
mv normal_sample_2.read1.fastq  RAWDATA/Sequencing data/Normal
mv normal_sample_2.read2.fastq  RAWDATA/Sequencing data/Normal
mv normal_sample_3.read1.fastq  RAWDATA/Sequencing data/Normal
mv normal_sample_3.read2.fastq  RAWDATA/Sequencing data/Normal
mv reference_chr22.fa Reference/
mv chr22.gtf Reference/

# Download Miniconda 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#install miniconda
~/miniconda3/bin/conda init bash
#create environment in conda
conda create --name bioinfo python=3.12.3
conda activate bioinfo

#Install tools
conda install -c bioconda <Toolname>
conda install -c bioconda fastqc
conda install -c bioconda Multiqc
conda install -c bioconda bwa 
conda install -c bioconda samtools


# Navigate to ReferenceFiles and index reference genome
cd RNAseq/Reference
bwa index reference_chr22.fa
samtools faidx reference_chr22.fa

# Use FastQC to assess the quality of each fastq file
cd /RNAseq/RAWDATA/Sequencing data/Normal
fastqc *.fastq
# Generate a MultiQC report to summarize the FastQC results
multiqc .

# and repeat for Cancer samples 

