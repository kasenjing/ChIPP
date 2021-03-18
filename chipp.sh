#!/bin/bash
start=$SECOND

print_help()
{
	echo ''
	echo -e "Usage: $0 -n [string] -f [DIR] -r [DIR] -q [INT] -p [broad or narrow]"
	echo -e "\t-n Data name you want to set"
	echo -e "\t-f Directory contain your raw fastq"
	echo -e "\t-r Directory contain bwa index"
	echo -e "\t-q Mapping quality"
	echo -e "\t-p Peak type"
}
# get arguments
while getopts n:f:r:q:p:h flag
do
	case "${flag}" in 
		n) name=${OPTARG};;
		f) fastq_dir=${OPTARG};;
		r) ref_dir=${OPTARG};;
		q) map_quality=${OPTARG};;
		p) peak_type=${OPTARG};;
		? | h) print_help
			exit 1;;
	esac
done


# make directory for processed data

if mkdir ${name}
then cd ${name}
else echo "directory has existed, change name"
	exit 1
fi

# collect the information of raw data
echo "calculating the # sequenceing"
fastq_seq_num=$(($(cat ${fastq_dir} | wc -l)/4 | bc))
echo "calculation is done."

#map reads to reference genome

echo "start mapping"
bwa aln ${ref_dir} ${fastq_dir} > ${name}.sai
bwa samse ${ref_dir} ${name}.sai ${fastq_dir} > ${name}.sam

# remove reads with quality smaller than -q specified INT

samtools view -h -q ${map_quality} ${name}.sam -F 07404 -o ${name}.q30.uni.sam

# calculate the number of seq do not pass quality control and -F 07404
header_num=$(($(awk '/^@/, /^[^@]/; /^[^@]/ {exit}' ${name}.q30.uni.sam | wc -l)-1 | bc))
total_num=$(($(cat ${name}.q30.uni.sam | wc -l)))
quality_align_num=$(($total_num-$header_num | bc))
echo "map completed"

# sort 

samtools sort -O BAM ${name}.q30.uni.sam -o ${name}.q30.uni.sorted.bam 

# mark duplicate reads

java -jar $Picard MarkDuplicates I=${name}.q30.uni.sorted.bam O=${name}.q30.uni.marked_dup.bam M=marked_dup_metrics.txt

# remove duplicated reads

samtools view -F 07404 ${name}.q30.uni.marked_dup.bam -o ${name}.q30.uni.dedup.bam

# calculate the num of duplicate

align_num_dedup=$(samtools view -c ${name}.q30.uni.dedup.bam)

# call peak with macs3

if [ ${peak_type} = "broad" ]
then
	macs3 callpeak -t ${name}.q30.uni.dedup.bam -g dm -n ${name} --broad
elif [ ${peak_type} = "narrow" ]
	macs3 callpeak -t ${name}.q30.uni.dedup.bam -g dm -n ${name} 
fi
echo -e "fastq_seq_num\t${fastq_seq_num}\nmapped_reads\t${mapped_reads}\nquality_align_num\t${quality_align_num}\nalign_num_dedup\t${align_num_dedup}\n" > statistic.txt
# make statistic plot
touch ${name}_stat.R
echo -e "counts <- c(${align_num_dedup},  ${quality_align_num}, ${fastq_seq_num})" >> ${name}_stat.R
echo -e "name <- c('deduplicated reads', 'quality_alignment','toral reads')" >> ${name}_stat.R
echo -e "barplot(counts, name.arg = name, horiz = TRUE, las = 1, space = c(0.1, 0.1, 0.1), width = c(0.15, 0.15, 0.15), main = 'reads statistic', axes = FALSE)" >> ${name}_stat.R
echo "duration=$((SECONDS - start))seconds"

