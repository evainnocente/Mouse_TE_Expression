!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=fastqc
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=result-%J.out
#SBATCH --error=error-%J.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL

## FastQC on RNAseq files

ml fastqc
mkdir trimmed_fastqc_output

for file in *.fastq.gz; do
  fastqc $file -o trimmed_fastqc_output/;
done

# End
