# workflow

This pipeline takes input files that are .fastq format(the file you usually get from sequencing company).

The first step of this pipeline is to map the sequences to reference genome using BWA.

While mapping reads to reference genome, there might be mismatches, deletions or insertions which cause low mapping quality.
Low quality means low reliablity. Also, one read can mapped to different positions in the genome because of the tolerance 
for mismatches, deletions or insertions. For example, read A is exactly the same with spot B while has one mismatch compared to spot C.
Mapping tool will retain the two alignments in the record file, so only one of the two record is needed for enrichment detection, obviously,
spot B is the best one. The file produced is .sam file.

The next step is to filter out the alignments with low quality, which is specify by users. Lower quality score will retain more alignments 
therefore more enrichment but cause low fidelity. Also, this step will retain one alignment for one read(spot B). When read A have several 
exactly same score(e.g. all perfectly match or all have one mismatch), one spot will be chosen randomly.

While preparing the sample for sequencing, multiple copies of one fragment may produced by PCR. They do not represent any of the enrichment of
target protein but are presented in the .sam file, they share be moved. Also, while sequencing a single amplification cluster, incorrectly 
detected as multiple clusters by the optical sensor of the sequencing instrument. They are called duplicates, the latter one is called optical 
duplicate.

The third step is to remove these duplicates. Alignments have exactly same start point and end point is regarded as duplicates. First, they are 
marked by using picard. Then they are moved with samtools.

The last step is to call peaks(enrichment) with MACS. This step produces file contain information about peak position and enrichment level.


