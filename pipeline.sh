#!/bin/bash
mkdir bioinformatics_assignment
cd bioinformatics_assignment
mkdir data
mkdir meta
mkdir logs
mkdir results

#download the fastq files and store them in untrimmed_fastq
cd data
mkdir untrimmed_fastq
cd untrimmed_fastq
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

#fastqc of untrimmed filescd ~/bioinformatics_assignment/data/untrimmed_fastq
cd ~/bioinformatics_assignment/data/untrimmed_fastq
fastqc NGS1.R1.fastq.gz
fastqc NGS1.R2.fastq.gz
mkdir ~/bioinformatics_assignment/results/fastqc_untrimmed_fastq
mv *fastqc* ~/bioinformatics_assignment/results/fastqc_untrimmed_fastq


#trimming using trimmomatic
mkdir $HOME/bioinformatics_assignment/data/trimmed_fastq
input1=$HOME/bioinformatics_assignment/data/untrimmed_fastq/NGS1.R1.fastq.gz
input2=$HOME/bioinformatics_assignment/data/untrimmed_fastq/NGS1.R2.fastq.gz
output1_paired="NGS1.R1.trimmed.paired.fq.gz"
output1_unpaired="NGS1.R1.trimmed.unpaired.fq.gz"
output2_paired="NGS1.R2.trimmed.paired.fq.gz"
output2_unpaired="NGS1.R2.trimmed.unpaired.fq.gz"

java -jar $HOME/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/trimmomatic.jar PE -phred33 -threads 4 \
"$input1" "$input2" \
"$output1_paired" "$output1_unpaired" \
"$output2_paired" "$output2_unpaired" \
ILLUMINACLIP:~/adapters/NexteraPE-PE.fa:2:30:10 \
TRAILING:25 MINLEN:50

mv ~/*fq.gz ~/bioinformatics_assignment/data/trimmed_fastq

#fastqc on the trimmed files
cd ~/bioinformatics_assignment/data/trimmed_fastq
fastqc NGS1.R1.trimmed.paired.fq.gz
fastqc NGS1.R2.trimmed.paired.fq.gz
cd

#alignment using BWA-mem
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mkdir ~/bioinformatics_assignment/data/reference
mv hg19.fa.gz ~/bioinformatics_assignment/data/reference
bwa index ~/bioinformatics_assignment/data/reference/hg19.fa.gz
mkdir ~/bioinformatics_assignment/data/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 \
~/bioinformatics_assignment/data/reference/hg19.fa.gz \
~/bioinformatics_assignment/results/trimmed_fastq/NGS1.R1.trimmed.paired.fq.gz \
~/bioinformatics_assignment/results/trimmed_fastq/NGS1.R2.trimmed.paired.fq.gz > \
~/bioinformatics_assignment/data/aligned_data/NGS1.hg19.sam

#convert to bam, sort and index with samtools
cd ~/bioinformatics_assignment/data/aligned_data
samtools view -h -b NGS1.hg19.sam > NGS1.hg19.bam
samtools sort NGS1.hg19.bam > NGS1.hg19.sorted.bam
samtools index NGS1.hg19.sorted.bam

#post alignment QC
#mark duplicates with picard
cd ~/bioinformatics_assignment/data/aligned_data
picard MarkDuplicates I=NGS1.hg19.sorted.bam O=NGS1.hg19.sorted.marked.bam M=marked_dup_metrics.txt
samtools index NGS1.hg19.sorted.marked.bam

#filter BAM based on mapping quality and bitwise flags
samtools view -F 1796 -q 20 -o NGS1.hg19.sorted.filtered.bam NGS1.hg19.sorted.marked.bam
samtools index NGS1.hg19.sorted.filtered.bam

#generate standard alignment statistics
samtools flagstat NGS1.hg19.sorted.filtered.bam
samtools idxstats NGS1.hg19.sorted.filtered.bam
samtools stats NGS1.hg19.sorted.filtered.bam


#variant calling
#convert the reference to text format
zcat ~/bioinformatics_assignment/data/reference/hg19.fa.gz > ~/bioinformatics_assignment/data/reference/hg19.fa

#index the reference
samtools faidx ~/bioinformatics_assignment/data/reference/hg19.fa

#call variants with freebayes
freebayes --bam ~/bioinformatics_assignment/data/aligned_data/NGS1.hg19.sorted.filtered.bam \
--fasta-reference ~/bioinformatics_assignment/data/reference/hg19.fa --vcf ~/bioinformatics_assignment/results/NGS1.hg19.vcf

#compress the vcf file
bgzip ~/bioinformatics_assignment/results/NGS1.hg19.vcf

#index the vcf
tabix -p vcf ~/bioinformatics_assignment/results/NGS1.hg19.vcf.gz

#filter the vcf
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/bioinformatics_assignment/results/NGS1.hg1S1.hg19.vcf.gz > ~/bioinformatics_assignment/results/NGS1.hg19.filtered.vcf

#intresect with bed file
bedtools intersect -header -wa -a ~/bioinformatics_assignment/results/NGS1.hg19.filtered.vcf \
-b ~/bioinformatics_assignment/data/untrimmed_fastq/annotation.bed > ~/bioinformatics_assignment/results/NGS1.hg19.filtered.intersect.vcf
bgzip ~/bioinformatics_assignment/results/NGS1.hg19.filtered.intersect.vcf
tabix -p vcf ~/bioinformatics_assignment/results/NGS1.hg19.filtered.intersect.vcf.gz

#annotation
#convert vcf to annovar input format
perl ~/annovar.latest/annovar/convert2annovar.pl \
-format vcf4 ~/bioinformatics_assignment/results/NGS1.hg19.filtered.intersect.vcf.gz > bioinformatics_assignment/results/NGS1.hg19.f>

#generate csv file
perl ~/annovar.latest/annovar/table_annovar.pl bioinformatics_assignment/results/NGS1.hg19.filtered.intersect.avinput \
~/annovar.latest/annovar/humandb/ -buildver hg19 -out ~/bioinformatics_assignment/results/NGS1.hg19.filtered.annovar \
-remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f \
-otherinfo -nastring . -csvout
