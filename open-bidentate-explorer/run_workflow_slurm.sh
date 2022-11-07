#!/bin/sh
#
#SBATCH --job-name="RUN_CREST"
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --account=education-as-msc-ce


python -m run_workflow.py > log_file_obe.out