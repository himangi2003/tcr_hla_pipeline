#!/bin/bash
#SBATCH --job-name=tcr_hla_pipeline
#SBATCH --output=logs/tcr_hla_%j.out
#SBATCH --ntasks=2
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --partition=campus-new

source ~/.bashrc
conda activate py311

python main.py
