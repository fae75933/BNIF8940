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

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi

module load BLAST+/2.9.0-gompi-2019b


blastn -num_threads 4 -query $DATADIR/sample.fasta -db /db/ncbiblast/nt/06042020/nt -out $OUTDIR/sample.fa.blastn.${SLURM_JOB_ID}.tsv -outfmt 6 -max_target_seqs=2
