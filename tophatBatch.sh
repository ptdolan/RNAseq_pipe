#!/bin/bash

for i in $SCRATCH/*R1*.fastq; do
echo $i
echo ${i/R1/R2}

cat <<EOM > ${i/R1_001.fastq/''}.slurm
#!/bin/bash
#SBATCH --job-name=tophatbatch
#
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G

set -o errexit

module purge
module load biology
module load bowtie2
module load tophat

tophat -o ${i/R1_001.fastq/_batchOutput} $SCRATCH/GRCm38/GRCm38_genome $i ${i/R1/R2}

EOM
done;
