#!/bin/bash
#SBATCH --job-name=testBLAST1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --output=/work/BNIF8940/fae75933/log.%j
#SBATCH --mail-user=fae75933@uga.edu
#SBATCH --mail-type=END,FAIL

OUTDIR="/work/gene8940/fae75933/"
DATADIR="/work/gene8940/instructor_data"
