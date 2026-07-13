#!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=telescope
#SBATCH --time=24:00:00 # increased
#SBATCH --mem=128G # increased
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL


BAMDIR="../08_Telescope_quantification/Sorted_bamfiles_output/"
ANNOTATIONS="../06_Index_with_RepeatMasker/repeatmasker_track/TE_repeatmasker_track.gtf" # Only TE annotations
OUTPUTDIR="./Telescope_output/"
LOGFILE="./logs/telescope_log.txt"

# Load necessary modules, required older versions
ml StdEnv/2020
ml python/3.7.9

# Activate the Python virtual environment that was made following instructions from https://github.com/finniej/telescope_noconda
source ../telescope_installation/telescope/telescope-env-py39/bin/activate

echo "$(date): virtual env activated" >> "$LOGFILE"

# Loop through all BAM files in bamdir
for file in "$BAMDIR"*.bam; do
    # Get base name without path or extension
    SAMPLENAME=$(basename "$file" _readnamesorted.bam)

    echo "Processing $SAMPLENAME" >> "$LOGFILE";
    # Need to make output directories first
    mkdir -p "$OUTPUTDIR"/"telescope_${SAMPLENAME}"

    # Run Telescope
    telescope assign "$file" "$ANNOTATIONS" \
    # Verbose output
    --debug \
    --logfile $LOGFILE \
    --outdir "$OUTPUTDIR/telescope_${SAMPLENAME}" \
    # Looks at this field in the GTF file to determine what the TE is
    --attribute transcript_id
    echo "$(date): $SAMPLENAME count completed" >> "$LOGFILE"
done

# Deactivate the environment (optional)
deactivate
