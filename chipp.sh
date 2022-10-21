#!/bin/bash
start=$SECOND

print_help()
{
	echo ''
	echo -e "Usage: $0 -n [string] -f [fastq DIR] -m [mate fastq DIR] -r [reference genome DIR] -q [INT] -p [broad or narrow]"
	echo -e "\t-n Data name you want to set"
	echo -e "\t-f Directory contain your raw fastq"
	echo -e "\t-r Directory contain bwa index"
	echo -e "\t-q Mapping quality"
	echo -e "\t-p Peak type"
	echo -e "\t-m mate fastq if not specified, run as single end mode."
}
# get arguments
while getopts n:f:r:q:p:m:h flag
do
	case "${flag}" in 
		n) name=${OPTARG};;
		f) fastq_dir=${OPTARG};;
		r) ref_dir=${OPTARG};;
		q) map_quality=${OPTARG};;
		p) peak_type=${OPTARG};;
		m) mate_fastq_dir=${OPTARG};;
		? | h) print_help
			exit 1;;
	esac
done


# make directory for processed data
pwd
cur=`pwd`
final_dir=${cur}/${name}
echo "${mate_fastq_dir}"

if mkdir $final_dir
then cd $final_dir
else echo "directory has existed, change name"
	exit 1
fi

# collect the information of raw data
echo "calculating the # sequenceing"
fastq_seq_num=$(($(zcat ${fastq_dir} | wc -l)/4 | bc))
echo "calculation is done."

#map reads to reference genome

echo "start mapping"
if [[ -z ${mate_fastq_dir} ]]
then
	echo "Detected single end input."
	bwa aln -t 16 ${ref_dir} ${fastq_dir} > ${name}.sai
	bwa samse ${ref_dir} ${name}.sai ${fastq_dir} > ${name}.sam
	# remove reads with quality smaller than -q specified INT.
	samtools view -h -q ${map_quality} ${name}.sam -F 1804 -o ${name}.q30.uni.sam
	samtools sort -O BAM ${name}.q30.uni.sam -o ${name}.q30.uni.sorted.bam
else
	if [[ -f ${mate_fastq_dir} ]]
	then
		echo "Detected pair end input."
		bwa mem -t 16 -M ${ref_dir} ${fastq_dir} ${mate_fastq_dir} > ${name}.sam
		# remove reads with quality smaller than -q specified INT in pair mode
		samtools view -h -q ${map_quality} -F 1804 -f 2 ${name}.sam -o ${name}.q30.uni.sam
		# remove orphan reads
		samtools fixmate -r ${name}.q30.uni.sam -o ${name}.q30.uni.tmp.sam
		samtools view -F 1804 -f 2 -u ${name}.q30.uni.tmp.sam | samtools sort -O BAM -o ${name}.q30.uni.sorted.bam	
	fi
fi

# calculate the number of seq do not pass quality control and -F 1804
quality_align_num=`samtools view -c ${name}.q30.uni.sorted.bam`
echo "map completed"


# mark duplicate reads

java -jar $picard MarkDuplicates I=${name}.q30.uni.sorted.bam O=${name}.q30.uni.dedup.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true

# remove duplicated reads

samtools view -F 07404 ${name}.q30.uni.marked_dup.bam -o ${name}.q30.uni.dedup.bam

# calculate the num of duplicate

align_num_dedup=$(samtools view -c ${name}.q30.uni.dedup.bam)

# call coverage with deeptools

samtools index ${name}.q30.uni.dedup.bam
bamCoverage -bs 1 -of bedgraph -b ${name}.q30.uni.dedup.bam -o ${name}.bedgraph
samtools flagstat ${name}.q30.uni.dedup.bam > ${name}.stat
# call peak with macs2

if [ ${peak_type} = "broad" ]
then
	echo -e "calling broad peak";
	macs2 callpeak -t ${name}.q30.uni.dedup.bam -n ${name} --broad
elif [ ${peak_type} = "narrow" ]
then
	echo -e "calling narrow peak";
	macs2 callpeak -t ${name}.q30.uni.dedup.bam -n ${name} 
fi

echo -e "fastq_seq_num\t${fastq_seq_num}\nquality_align_num\t${quality_align_num}\nalign_num_dedup\t${align_num_dedup}\n" > statistic.txt
# make statistic plot
#touch ${name}_stat.R
#echo -e "counts <- c(${align_num_dedup},  ${quality_align_num}, ${fastq_seq_num})" >> ${name}_stat.R
#echo -e "name <- c('deduplicated reads', 'quality_alignment','toral reads')" >> ${name}_stat.R
#echo -e "barplot(counts, name.arg = name, horiz = TRUE, las = 1, space = c(0.1, 0.1, 0.1), width = c(0.15, 0.15, 0.15), main = 'reads statistic', axes = FALSE)" >> ${name}_stat.R
#echo "duration=$((SECONDS - start))seconds"

