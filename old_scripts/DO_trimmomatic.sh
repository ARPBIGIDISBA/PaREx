SAMPLES=$1
LINEAGE=$2

mkdir /home/micro/Analysis/Trimmomatic/"$LINEAGE" 

while read line
do

	mkdir /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"  


		java -jar /home/micro/Programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /home/micro/Analysis/assemblies/fastq-MiSeq/$LINEAGE/"$line"_R1_001.fastq.gz /home/micro/Analysis/assemblies/fastq-MiSeq/$LINEAGE/"$line"_R2_001.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.1P.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.1U.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.2P.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.2U.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 


	mv /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.1P.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"_trim_R1_001.fastq.gz
	
	mv /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"/"$line".trimmed.2P.fastq.gz /home/micro/Analysis/Trimmomatic/"$LINEAGE"/"$line"_trim_R2_001.fastq.gz


	
done < $SAMPLES



	
