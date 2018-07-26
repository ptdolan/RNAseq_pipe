 This directory contains the scripts used to map and analyze RNAseq data in batch on a cluster (specifically written with for slurm job submission on sherlock.stanford.edu. 

Instructions:

---------
A. TOPHAT
---------

1. Move all fastq files to $SCRATCH on cluster
2. From within $SCRATCH directory (or dir where fastq's live), run tophatBatch.sh.
  > sh /path/to/tophatBatch.sh

This will generate all of the slurm command files for the directory listed in the *Batch.sh file.

3. Then run:
  > for i in *tophat.slurm; do
  > sbatch $i
  > done

------------
B. CUFFLINKS
------------

1. Confirm that all tophat2 runs have completed, 'align_summary.txt' and 'mapped_reads.bam' should now be in the "*batchOutput/" directories where * is the sample ID info. 

2. From within $SCRATCH (where fastq's live), run Batch.sh.
> sh /path/to/cufflinksBatch.sh

This will generate all of the slurm command files for cufflinks based on the folders from tophat.

3. Then run:

> for i in *cufflinks.slurm; do
> sbatch $i
> done





