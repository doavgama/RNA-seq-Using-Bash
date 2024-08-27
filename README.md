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
and repeat for Cancer samples 

#Error trimming to enhance the quality 
#overrepresented sequence trimming
install BBMap
conda install -c bioconda BBMap

# to remove over represented seq
bbduk.sh in1=normal_sample_1.read1.fastq in2=normal_sample_1.read2.fastq \
         out1=trimmed_R1.fastq out2=trimmed_R2.fastq \
         ref=read2.fasta k=31 hdist=1
bbduk.sh in1=normal_sample_2.read1.fastq in2=normal_sample_2.read2.fastq \
         out1=trimmednormals2r1.fastq out2=trimmednormals2r2.fastq \
         ref=read2.fasta k=31 hdist=1
bbduk.sh in1=normal_sample_3.read1.fastq in2=normal_sample_3.read2.fastq \
         out1=trimmednormals3r1.fastq out2=trimmednormals3r2.fastq \
         ref=read2.fasta k=31 hdist=1
 ref:FASTA file containing overrepresented sequences to be removed.
 
 PS : i had to gather all the overrepresented sequence manually from the fastqc report  
 
 check the quality again 
 fastqc trimmed_R1.fastq trimmed_R2.fastq trimmednormals2r1 trimmednormals2r2 trimmednormals3r1.fastq trimmednormals3r2.fastq

#  trim specifc seq 
seqtk trimfq -b 18 trimmed_s1R1.fastq > trimmnormals1r1.fastq
seqtk trimfq -b 18 trimmed_s1R2.fastq > trimmnormals1r2.fastq


seqtk trimfq -b 18 trimmednormals2r1.fastq > trimms2r1.fastq
seqtk trimfq -b 18 trimmednormals2r2.fastq > trimms2r2.fastq

seqtk trimfq -b 18 trimmednormals3r1.fastq > trimms3r1.fastq
seqtk trimfq -b 18 trimmednormals3r2.fastq > trimms3r2.fastq

check quality 
fastqc trimms1r1.fastq trimms1r2.fastq trimms2r1.fastq trimms2r2.fastq trimms3r1.fastq trimms3r2.fastq


# Remove overrepresented seq from cancer samples 

 bbduk.sh in1=cancer_sample_1.read1.fastq in2=cancer_sample_1.read2.fastq          out1=trimmed_cancers1r1.fastq out2=trimmed_cancers1r2.fastq          ref=over1.fasta k=31 hdist=1

bbduk.sh in1=cancer_sample_2.read1.fastq in2=cancer_sample_2.read2.fastq          out1=trimmed_cancers2r1.fastq out2=trimmed_cancers2r2.fastq          ref=over2.fasta k=31 hdist=1

bbduk.sh in1=cancer_sample_3.read1.fastq in2=cancer_sample_3.read2.fastq          out1=trimmed_cancers3r1.fastq out2=trimmed_cancers3r2.fastq          ref=over2.fasta k=31 hdist=1

# Trim the first 18 seq

seqtk trimfq -b 18 trimmed_cancers1r1.fastq > trimms1r1.fastq
seqtk trimfq -b 18 trimmed_cancers1r2.fastq > trimms1r2.fastq

seqtk trimfq -b 18 trimmed_cancers2r1.fastq > trimms2r1.fastq
seqtk trimfq -b 18 trimmed_cancers2r2.fastq > trimms2r2.fastq

seqtk trimfq -b 18 trimmed_cancers3r1.fastq > trimms3r1.fastq
seqtk trimfq -b 18 trimmed_cancers3r2.fastq > trimms3r2.fastq


# check quality 

fastqc trimms1r1.fastq trimms1r2.fastq trimms2r1.fastq trimms2r2.fastq trimms3r1.fastq  trimms3r2.fastq

#Move all files to Reference directory 
mv trimms1r1.fastq trimms1r2.fastq trimms2r1.fastq trimms2r2.fastq trimms3r1.fastq trimms3r2.fastq RNAseq/References/cancer/
mv trimms1r1.fastq trimms1r2.fastq trimms2r1.fastq trimms2r2.fastq trimms3r1.fastq trimms3r2.fastq RNAseq/References/normal/

# Reference Genome Preparation
o Index the reference genome using BWA
# install the tool 
conda install -c bioconda BWA
# indexing 
cd RNAseq/References
bwa index reference_chr22.fa


# Align the reads to the reference genome using BWA MEM
cd RNAseq/References/normal
bwa mem -t 4 reference_gemome_chr22.fa trimmnormals1r1.fastq trimmnormals1r2.fastq > normaloutput1.sam
bwa mem -t 4 reference_gemome_chr22.fa trimmnormals2r1.fastq trimmnormals2r2.fastq > normaloutput2.sam
bwa mem -t 4 reference_gemome_chr22.fa trimmnormals3r1.fastq trimmnormals3r2.fastq > normaloutput3.sam
# cancer 
cd - (to get back to the last dir)
cd cancer
bwa mem -t 4 reference_gemome_chr22.fa trimms1r1.fastq trimms1r2.fastq > canceroutput1.sam
bwa mem -t 4 reference_gemome_chr22.fa trimms2r1.fastq trimms2r2.fastq > canceroutput2.sam
bwa mem -t 4 reference_gemome_chr22.fa trimms3r1.fastq trimms3r2.fastq > canceroutput3.sam

# Convert the resulting SAM files to sorted BAM files using Samtools
samtools view -S -b normaloutput1.sam | samtools sort -o sorted_normaloutput1.bam
samtools view -S -b normaloutput2.sam | samtools sort -o sorted_normaloutput2.bam
samtools view -S -b normaloutput3.sam | samtools sort -o sorted_normaloutput3.bam

samtools view -S -b canceroutput1.sam | samtools sort -o sorted_canceroutput1.bam
samtools view -S -b canceroutput2.sam | samtools sort -o sorted_canceroutput2.bam
samtools view -S -b canceroutput3.sam | samtools sort -o sorted_canceroutput3.bam

# Index the sorted BAM files
samtools index sorted_normaloutput1.bam
samtools index sorted_normaloutput2.bam
samtools index sorted_normaloutput3.bam

samtools index sorted_canceroutput1.bam
samtools index sorted_canceroutput2.bam
samtools index sorted_canceroutput3.bam

#Read Quantification
#Use featureCounts to quantify reads mapped to genes

install featureCounts 
sudo apt subread
# Generate a count matrix for all samples

featureCounts -T 4 -a chr22.gtf -o normalcounts.txt -p sorted_normaloutput1.bam sorted_normaloutput2.bam sorted_normaloutput3.bam
featureCounts -T 4 -a chr22.gtf -o cancercounts.txt -p sorted_canceroutput1.bam sorted_canceroutput2.bam sorted_canceroutput3.bam

PS : i had to re-trim all my files to fix this error : No paired-end reads were detected in paired-end read library : sorted_(all my bam files either it's normal or cancer).
mv cancercounts.txt cancercounts.txt.summary RNAseq/Results
mv normalcounts.txt normalcounts.txt.summary RNAseq/Results

#Preview the Data
less cancercounts.txt.summary
less normalcounts.txt.summary 
head cancercounts.txt 
head normalcounts.txt
