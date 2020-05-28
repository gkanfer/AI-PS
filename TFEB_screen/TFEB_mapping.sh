#!/bin/sh

#SBATCH --job-name="Gil_41"
#SBATCH --partition=norm
#SBATCH --cpus-per-task=20
#SBATCH --time=50:00:00
#SBATCH --mem=100g
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:256
#SBATCH --job-name="Gil_41"       

module load bcl2fastq
module load umitools
module load trimmomatic
module load fastxtoolkit
module load bowtie/1.1.2
module load samtools
module load umitools


input_path=$(pwd)

mkdir tmp
bcl2fastq --no-lane-splitting --runfolder-dir ./190917_NB552201_0006_AHMMNVAFXY/ --output-dir tmp &> tmp/bcl2fastq.log 
bcl2fastq --no-lane-splitting --runfolder-dir $input_path --output-dir tmp &> tmp/bcl2fastq.log

cd tmp 

gunzip Undetermined_S0_R1_001.fastq.gz

umi_tools extract --stdin=Undetermined_S0_R1_001.fastq --extract-method='regex' --bc-pattern="(?P<umi_1>.{3}T{1}.{3}T{1}.{3}T{1}.{8})(?P<discard_1>GCACAAAAGGAAAC.*AGTATCCCTTGGAGAACCACCTTGTTGG){s<=2}(.*)(?P<discard_2>GTTTAAGAGCTAAGCT.*){s<=2}" --log=processed.log --stdout processed.fastq.gz 
java -jar $TRIMMOJAR SE -phred33 processed.fastq.gz trimmed_processed.fastq.gz SLIDINGWINDOW:4:15 MINLEN:19
gunzip trimmed_processed.fastq.gz
#fastq-sample -n 69000000 trimmed_processed.fastq -o trimmed_processed_sample
fastx_trimmer -f 1 -l 19  -i trimmed_processed.fastq -o trimmed_processed_fastx.fastq
#gzip trimmed_processed_fastx.fastq

path="/data/kanferg/DEEP_SEQ/refarance_extantion/Crispri_whoule_genome/Crispri_h"

#gunzip trimmed_processed_fastx.fastq.gz

bowtie --threads 10 --tryhard -n 0 -k 1 -l 19 $path -q trimmed_processed_fastx.fastq --chunkmbs 200 --sam | samtools view -Sb - > Post_Map_e.bam

samtools sort Post_Map_e.bam -o Sorted_Post_Map_e.bam
samtools index Sorted_Post_Map_e.bam

umi_tools count --per-contig -I Sorted_Post_Map_e.bam > e_UMI_countList_test_e.txt
