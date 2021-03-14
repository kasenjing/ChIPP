# ChIP-seq
This pipeline is designed for ChIP-seq data analysis. 
# Workflow

# dependencies

#### tools used in this pipeline
* bwa -version 0.7.17
 
        https://github.com/lh3/bwa
* samtools -version 1.5
 
        https://github.com/samtools/samtools
* picard

        https://github.com/broadinstitute/picard

#### add tools to path


	vi ~/.bashrc
	export PATH=$PATH:/home/kasenjing/biotools/bwa
	export Picard=/home/kasenjing/biotools/picard/build/libs/picard.jar
	source ~/.bashrc

# usage

	./chipp.sh -n [string] -f [DIR] -r [DIR] -q [INT] -p [broad or narrow]
	
Arguments | Description
----------|-----------
   -n | Name of the factor you are dueling with.
   -f | Directory contain the .fastq file.
   -r | Directory contain the bwa index file.
   -q | Filter out alignments with MAPQ < [INT].
   -p | Peak type.

