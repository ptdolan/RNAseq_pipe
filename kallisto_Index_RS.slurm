#!/bin/bash
#SBATCH --job-name=KallistoTranscriptIndex
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH -e KallistoIndexError_REFSEQ.out
module load biology
module load kallisto
which kallisto
 
kallisto index --make-unique -i $SCRATCH/GRCm38/GRCm38_refseq_transcripts.idx $SCRATCH/GRCm38/GCF_000001635.26_GRCm38.p6_rna.fna.gz

