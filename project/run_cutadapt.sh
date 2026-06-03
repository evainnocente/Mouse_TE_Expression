# Running cutadapt to trim and filter reads
# Create required virtual environment locally
ml python/3.11
ENVDIR=/tmp/$RANDOM
virtualenv --no-download $ENVDIR
source $ENVDIR/bin/activate
pip install --no-index --upgrade pip
pip install --no-index cutadapt # Nothing else required
pip freeze --local > cutadapt_requirements.txt # freezes versions in cutadapt file
cat cutadapt_requirements.txt # check, looks good
deactivate
rm -rf $ENVDIR # remove, success

# parameters
# -q 25
# -m 40 -length
# -o trimmed_$file -p trimmed_$file
# -u 6 -U 6 #positive number = from 5' end, trimming number of bases
# -O/--overlap 5, trim bases overlapping with adapters, stringency equivalent,  nucleotide overlap with Illumina primer sequence for trimming

# Ran following as a script submitted to compute nodes
#!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=cutadapt
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1 # specify 1 node to make things easier
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL

module load python/3.11 # default
virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r cutadapt_requirements.txt

# from within pwd
cutadapt --version
for read1 in *_R1_*.fastq.gz; do
  read2=${read1/_R1_/_R2_} # replacing R1 with R2 to iterate over both reads
  cutadapt -q 25 -m 40 -u 6 -U 6 --overlap 5 -o trimmed_${read1} -p trimmed_${read2} ${read1} ${read2}
done

deactivate

# End
