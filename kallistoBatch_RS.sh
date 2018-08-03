#!/bin/bash

for i in $SCRATCH/*R1*.fastq; do
echo $i
echo ${i/R1/R2}

cat <<EOM > ${i/R1_001.fastq/Kallisto_RS}.slurm
#!/bin/bash
#SBATCH --job-name=kallistoBatch
#
#SBATCH --time=20:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

set -o errexit

module purge
module load biology
module load kallisto
 
kallisto quant -t 4 -b 100 -o ${i/R1_001.fastq/_Kallisto_RS_Output} -i $SCRATCH/GRCm38/GRCm38_refseq_transcripts.idx $i ${i/R1/R2}

EOM
done;
