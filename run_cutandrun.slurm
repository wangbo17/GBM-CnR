#!/bin/bash
#SBATCH --job-name=cutandrun
#SBATCH --output=cutandrun.out
#SBATCH --error=cutandrun.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G

conda activate nf_env
nextflow run cutandrun.nf -resume 2eca8adf12