#!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=starannotation
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1 # specify 1 node to make things easier
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL

ml StdEnv/2023
ml star/2.7.11b

FASTADIR="./Ref_genome_fastas/Mus_musculus.GRCm39.dna.primary_assembly.fa" # Only 1 fasta now, not an array
OUTPUTDIR="./Genome_index_dir/"
GTFFILE="./GTF_file_for_indexing/combined_refgenome_TEs.gtf"
LOGFILE="./logfile/index_creationlog.txt"

echo "$(date): Directory of ref genome fastas = $FASTADIR, directory of ref genome annotation incl TEs = $GTFFILE, directory for output files = $OUTPUTDIR" >> "$LOGFILE"
echo "$(date): Starting index creation" >> "$LOGFILE"
STAR --runMode genomeGenerate --genomeDir "$OUTPUTDIR" --genomeFastaFiles "$FASTADIR" --sjdbGTFfile "$GTFFILE" --sjdbOverhang 149 #1 less than read length, removed --sjdbGTFtagExonParentTranscript Parent
echo "$(date): Genome index with transposable element repeat masker file and reference gtf including scaffolds complete" >> "$LOGFILE"

# End
