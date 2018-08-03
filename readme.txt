This directory contains the scripts used to map and analyze RNAseq data in batch on a cluster (specifically written with for slurm job submission on sherlock.stanford.edu. 
 
Two options currently, tophat2-bowtie2-cufflinks pipeline (slow), or the Kallisto-Sleuth pipeline (fast).

Instructions:

---------------------------------------
1A. BOWTIE INDICES
---------------------------------------

1. Bowtie requires an indexed representation of the target DNA to align to. To generate the index (here generating Long indices, see Bowtie Documentation for more info), edit the 'bowtieIndexL.sh' file to point to your target fasta/fa.gz.

2. Run indexing on slurm cluster:
	> sbatch bowtieIndexL.sh

---------------------------------------
1A. TOPHAT 
---------------------------------------

1. Move all fastq files to $SCRATCH on cluster

2. Edit tophatBatch.sh to point toward your input and outputs (index file, target fasta, fasts).

3. From within $SCRATCH directory (or dir where fastq's live), run tophatBatch.sh.
  > sh /path/to/tophatBatch.sh

This will generate all of the slurm command files for the directory listed in the *Batch.sh file.

4. Then run:
  > for i in *tophat.slurm; do
  > sbatch $i
  > done

---------------------------------------
1C. CUFFLINKS
---------------------------------------

1. Confirm that all tophat2 runs have completed, 'align_summary.txt' and 'mapped_reads.bam' should now be in the "*batchOutput/" directories where * is the sample ID info. 

2. From within $SCRATCH (where fastq's live), run Batch.sh.
  > sh /path/to/cufflinksBatch.sh

This will generate all of the slurm command files for cufflinks based on the folders from tophat.

3. Then run:

  > for i in *cufflinks.slurm; do
  > sbatch $i
  > done


#######################################

---------------------------------------
2A. Kallisto index - alignment and counting on Cluster
---------------------------------------

1. Kallisto requires an indexed representation of the target DNA to align to. To generate the index, edit the 'kallistoIndex.sh' file to point to your target fasta/fa.gz.

2. Run indexing on slurm cluster, like sherlock.stanford.edu:
  > sbatch kallistoIndex.sh


---------------------------------------
2B. Kallisto quant - alignment and counting on Cluster
---------------------------------------


Kallisto Manual: https://pachterlab.github.io/kallisto

1. Move all fastq files to $SCRATCH on cluster (faster I/O)
2. From within $SCRATCH directory (or dir where fastq's live), run 'tophatBatch.sh' .
  > sh /path/to/kallistoBatch.sh

This will generate all of the slurm command files for the directory listed in the *Batch.sh file.

3. Then run:
  > for i in *kallisto.slurm; do
  > sbatch $i
  > done

  
---------------------------------------
2C. Sleuth - analysis
---------------------------------------
Sleuth Manual: https://pachterlab.github.io/sleuth/about

1. Place all "KallistoOutput/" directories into same parent directory.
2. Edit SleuthAnalysis.R to match your file directory and experimental model structure. 
3. Run Sleuth to generate outputs.
4. To view interactive report of results (requires 'Shiny' package) type sleuth_live(SO), where 'SO' is your sleuth object. 



