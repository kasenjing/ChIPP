# ChIP-seq
This pipeline is designed for ChIP-seq data analysis. This pipeline is designed for short single-end reads, provides r scripts of reads number passes through each step for statistics.
# Workflow
[workflow](https://github.com/kasenjing/ChIPP/blob/main/CHIPPworkflow.pdf)

# dependencies

#### tools used in this pipeline
* bwa -version 0.7.17
 
	https://github.com/lh3/bwa
* samtools -version 1.5
 
	https://github.com/samtools/samtools
* picard

 	https://github.com/broadinstitute/picard
* macs

	https://github.com/macs3-project/MACS
	

#### add tools to path


	vi ~/.bashrc
	export PATH=$PATH:/home/kasenjing/biotools/	# directory contain the tools
	export Picard=/home/kasenjing/biotools/picard/build/libs/picard.jar
	source ~/.bashrc

# usage

	./chipp.sh -n [string] -f [DIR] -r [DIR] -q [INT] -p [broad or narrow]
	
Arguments | Description
----------|-----------
   -n | Name of the factor you are dealing with.
   -f | Directory contain the .fastq file.
   -r | Directory contain the bwa index file.
   -q | Filter out alignments with MAPQ < [INT].
   -p | Peak type.

