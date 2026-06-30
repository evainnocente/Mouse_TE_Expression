#!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=starannotation
#SBATCH --time=24:00:00 # increased
#SBATCH --mem=128G # increased
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL

ml star/2.7.11b
ml samtools/1.22.1

INDEXDIR="../06_Index_with_RepeatMasker/Genome_index_dir/" #same index
READSDIR="../Untrimmed/04_5_Alignment_rawreads/Raw_fastqs/" # assuming all reads are in this directory ## nchange
OUTPUTDIR="./Alignment_output/" # made beforehand
LOG_FILE="./logfiles/mapping_logfile.txt" # made beforehand
THREADS=8

echo "$(date): Starting STAR alignment on paired-end trimmed data" >> "$LOG_FILE"

for R1 in $READSDIR*_R1_*.fastq.gz; do  # loop through fwd reads
  R2="${R1/_R1_/_R2_}"; # Create variable of reverse reads
  SAMPLENOPREFIX="${R1/trimmed_/}"; # Remove prefix
  SAMPLENOR1="${SAMPLENOPREFIX/_R1/}"; # Remove R1
  SAMPLENAME=$(basename $SAMPLENOPREFIX .fastq.gz); # Create sample name by removing extension

  echo "$(date): Aligning sample ${SAMPLENAME}" >> "$LOG_FILE";

  STAR --runThreadN $THREADS \
    --genomeDir "$INDEXDIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUTPUTDIR}/${SAMPLENAME}_untrimmed_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within KeepPairs \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 200 \
    --outFilterMatchNminOverLread 0.3 \
    --outFilterMatchNmin 30 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignSJoverhangMin 5 \
    --alignSJDBoverhangMin 3 \
    --alignIntronMin 20 \
    --sjdbScore 2 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --outSAMattributes Standard \
    --outFilterType BySJout

# module avail samtools
# samtools/1.22.1 (bio,D)

# Index the output BAM file
  BAM_FILE=${OUTPUTDIR}${SAMPLENAME}_untrimmed_Aligned.sortedByCoord.out.bam
  echo "$(date): Indexing BAM file ${BAM_FILE}" >> "$LOG_FILE"
  samtools index "$BAM_FILE"

# Log completion of sample
  echo "$(date): Finished alignment and indexing of ${SAMPLENAME}" >> "$LOG_FILE"

done
