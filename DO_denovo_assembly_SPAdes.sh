SAMPLES=$1
LINEAGE=$2
PATH =/home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/ 

mkdir $PATH/$LINEAGE

while read line
do

	mkdir $PATH/$LINEAGE/$line

	python3 /home/micro/Programs/SPA3/spades-3.15.5/assembler/spades.py -o $PATH/$LINEAGE/$line -1 /home/micro/Analysis/assemblies/fastq-Miseq/$LINEAGE/"$line"_L001_R1_001.fastq -2 /home/micro/Analysis/assemblies/fastq-Miseq/$LINEAGE/"$line"_L001_R2_001.fastq --careful

	mv $PATH/$LINEAGE/"$line"/contigs.fasta $PATH/$LINEAGE/"$line".SPAdes.denovoassembly.fasta
	
done <$SAMPLES
