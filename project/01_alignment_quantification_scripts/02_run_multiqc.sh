# Script to run MultiQC on FastQC output

# First load modules locally to make text list of software that the python environment will require
python/3.11
ENVDIR=/tmp/$RANDOM
virtualenv --no-download $ENVDIR
source $ENVDIR/bin/activate
pip install --no-index --upgrade pip
ml gcc arrow/24.0.0 # Required for this install
pip install --no-index multiqc # Check version
# Capture env requirements
pip freeze --local > multiqc_requirements.txt # Freezes versions
cat multiqc_requirements.txt # check, looks good
deactivate
rm -rf $ENVDIR # Remove environment locally

# Start script

!/bin/bash
#SBATCH --account=def-jlamarre
#SBATCH --job-name=multiqc
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1 # specify 1 node to make things easier
#SBATCH --output=result-%J.out
#SBATCH --error=error-%J.err
#SBATCH --mail-user=einnocen@uoguelph.ca
#SBATCH --mail-type=ALL

module load python/3.11 # default
virtualenv --no-download $SLURM_TMPDIR/env # Activate environment within job as recommended
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r multiqc_requirements.txt # Use pre-existing list of requirements for install

# From within pwd
multiqc --version
multiqc -v . -o trimmed_multiqc_output/ # Verbose
deactivate

# End script
