SAMPLES=$1
LINEAGE=$2

mkdir /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/$LINEAGE 

while read line
do

	mkdir /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/$LINEAGE/$line  


	python /home/micro/Programs/SPAdes-3.15.0/bin/spades.py -o /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/$LINEAGE/$line -1 /home/micro/Analysis/assemblies/fastq-MiSeq/$LINEAGE/"$line"_R1_001.fastq -2 /home/micro/Analysis/assemblies/fastq-MiSeq/$LINEAGE/"$line"_R2_001.fastq --careful


	mv /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/$LINEAGE/$line/contigs.fasta /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/$LINEAGE/"$line".SPAdes.denovoassembly.fasta

		
done < $SAMPLES


