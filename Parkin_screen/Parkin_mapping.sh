#!/bin/sh

#SBATCH --job-name="Park"
#SBATCH --partition=norm
#SBATCH --cpus-per-task=20
#SBATCH --time=05:00:00
#SBATCH --mem=100g
#SBATCH --mail-type=ALL
#SBATCH --gres=scratch:256
#SBATCH --job-name="Park"


module load bcl2fastq
module load umitools
module load trimmomatic
module load fastxtoolkit
module load bowtie/1.1.2
module load samtools
module load umitools
module load bedtools/2.25.0

input_folder=$pwd
folder=$1
file=$2
ref="./refarance_extantion/Crispri_h1"

mkdir $folder
cp $file ./$folder

cd $folder

mv $file processed.fastq

gzip processed.fastq

java -jar $TRIMMOJAR SE -phred33 processed.fastq.gz trimmed_processed.fastq.gz SLIDINGWINDOW:4:15 MINLEN:19
gunzip trimmed_processed.fastq.gz
fastx_trimmer -f 2 -l 19  -i trimmed_processed.fastq -o trimmed_processed_fastx.fastq
bowtie --threads 10 --tryhard -n 0 -k 1 -l 19 $ref -q trimmed_processed_fastx.fastq --chunkmbs 200 --sam | samtools view -Sb - > f_trim.bam
bedtools bamtobed -i f_trim.bam > f_trim.bed
cut -f 1 f_trim.bed|sort|uniq -c > ind1.bed
