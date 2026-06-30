#!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=adaptertrim
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL


ml fastp/1.0.1 # default

READSDIR="./Untrimmed_fastqs/" # assuming all reads are in this directory
OUTPUTDIR="./Adapters_trimmed_reads_output/" # made beforehand
LOG_FILE="./logfiles/trimming_log.txt" # made beforehand

echo "$(date): Starting trimming on paired-end data" >> "$LOG_FILE"

for read1 in $READSDIR*_R1_*.fastq.gz; do
  read2=${read1/_R1_/_R2_}; # replacing R1 with R2 to iterate over both reads
  echo "$(date) Finished trimming on $read1 and $read2" >> "$LOG_FILE";
  fastp -q 25 -l 40 --trim_front1 6 --trim_front2 6 -i $read1 -o ${OUTPUTDIR}adapter_trimmed_$(basename {read1}) -I $read2 -O ${OUTPUTDIR}adapter_trimmed_$(basename {read2}); #trims only bases on 5' ends of each read
  echo "$(date) Finished trimming on $read1 and $read2" >> "$LOG_FILE";
done

echo "$(date): Finished trimming on paired-end data" >> "$LOG_FILE"

# End
