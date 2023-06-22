#!/bin/bash

SAMPLES=$1 #load a list of files 

mkdir /home/micro/Analysis/assemblies/fastq-Miseq/CNB
mkdir /home/micro/Analysis/sam_files/CNB
#output de salida donde se generarán los sam comprimidos /media/micro/easystore/sam_files/mutantes/

while read line
do

        echo "Aligning $line..."
#1
#cp /media/micro/2F1DC2EC2D57481B/CNB_CAZ_CAZAVI/"$line"_1.fastq.gz /home/micro/Escritorio/prova/
#cp /media/micro/Expansion/prova/"$line"_2.fastq.gz /home/micro/Escritorio/prova/

#2
#gzip -d /home/micro/Escritorio/prova/"$line"_1.fastq.gz
#gzip -d /home/micro/Escritorio/prova/"$line"_2.fastq.gz


#1
cp /media/micro/2F1DC2EC2D57481B/CNB_CAZ_CAZAVI/"$line"_1.fastq /home/micro/Analysis/assemblies/fastq-Miseq/CNB
cp /media/micro/2F1DC2EC2D57481B/CNB_CAZ_CAZAVI/"$line"_2.fastq /home/micro/Analysis/assemblies/fastq-Miseq/CNB

#3
	/home/micro/Programs/bowtie2-2.2.6/bowtie2 --phred33 -x /home/micro/indexed_libs/PAO1/PAO1 -q -1 /home/micro/Analysis/assemblies/fastq-Miseq/CNB/"$line"_1.fastq -2 /home/micro/Analysis/assemblies/fastq-Miseq/CNB/"$line"_2.fastq -X 1000 -S /home/micro/Analysis/sam_files/CNB/"$line"_mapPAO1.sam

#4

rm /home/micro/Analysis/assemblies/fastq-Miseq/CNB/"$line"_1.fastq
rm /home/micro/Analysis/assemblies/fastq-Miseq/CNB/"$line"_2.fastq

#5

gzip /home/micro/Analysis/sam_files/CNB/"$line"_mapPAO1.sam 

#6

mv /home/micro/Analysis/sam_files/CNB/"$line"_mapPAO1.sam.gz /media/micro/easystore/sam_files/mutantes/

echo "Done aligning $line..."

done < $SAMPLES

#3
	/home/micro/Programs/bowtie2-2.2.6/bowtie2 --phred