#!/bin/bash
#SBATCH --job-name=homework1
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --output=/work/BNIF8940/fae75933/log.%j
#SBATCH --mail-user=fae75933@uga.edu
#SBATCH --mail-type=END,FAIL

OUTDIR="/work/gene8940/fae75933/"
cd $OUTDIR
curl -s ftp://ftp.ensemblgenomes.org/pub/bacteria/release-37/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz | gunzip -c > $OUTDIR/ecoli_MG1655.gff
grep -c "CDS" ecoli_MG1655
