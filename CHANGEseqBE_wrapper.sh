#!/bin/bash
#SBATCH --job-name IGU_CHANGEseq
#SBATCH -n 32
#SBATCH -o %j.out
#SBATCH -e %j.err


## ACTIVATE CONDA
eval "$(conda shell.bash hook)"
conda activate changeseqbe
SCRIPT_PATH=/home/thudson/projects/changeseq_be/CHANGEseqBE_wrapper.sh
# Get the directory where the script is located
SCRIPT_DIR="$(dirname "$(realpath "$SCRIPT_PATH")")"

# Now you can reference this directory in your script
#echo "parameters manifest(csv): $1"
#echo "samples manifest(csv): $2"
echo "manifest(csv): $1"
echo "command: $2"
echo "sample: $3" # all or a specific sample
echo "making yaml file"

cd $SCRIPT_DIR/changeseq
python $SCRIPT_DIR/changeseq/changeseq.py $2 --manifest $1 --sample $3


