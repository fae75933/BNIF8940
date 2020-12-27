#!/bin/bash
#SBATCH --job-name=optimization
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH --output=/work/gene8940/fae75933/log.%j
#SBATCH --mail-user=fae75933@uga.edu
#SBATCH --mail-type=END,FAIL

#path to files: cd /scratch/zlewis/Run119/FastQs/Ncrassa


OUTDIR="/work/gene8940/fae75933/optimization"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

Threads=6
MEMORY=24

FILES="/scratch/zlewis/Run119/FastQs/Ncrassa/*R1_001.fastq.gz" #Don't forget the *

OUTDIR="/scratch/fae75933/CUT&RUN"

mkdir "$OUTDIR/SortedBamFiles"
mkdir "$OUTDIR/Bigwigs"
mkdir "$OUTDIR/MACSout"

module load deepTools/3.3.1-intel-2019b-Python-3.7.4
ml BWA/0.7.17-GCC-8.3.0
# ml SAMtools/1.9-foss-2016b
ml SAMtools/1.9-GCC-8.3.0

#Iterate over the files in the fastq folder and perform desired analysis steps
for f in $FILES
do


##define the variable $file to extract just the filename from each input file. Note that the variable $f will contain the entire path. Here you will extract the name of the file from the path and use this to name files that descend from this original input file.

file=${f##*/}

read2=$(echo "$f" | sed 's/R1.fastq.gz/R2.fastq.gz/g') #this is a variable that stores the name of the read 2 file for matching paired-end fastq files. !!!!!!if this doesn't work it is probably because the end of the fastq file name has a different format. Make sure it matches.

#filename variable for the sorted bam file
sorted="$OUTDIR/SortedBamFiles/$file"

#filename variable for the deeptools big wig output
bigwig="$OUTDIR/Bigwigs/$file.bw"

bwa mem -M -v 3 -t $THREADS /scratch/fae75933/genomesfolder/GCA_000182925_neurospora.fna $f $read2 | samtools view -bhu - | samtools sort -T $file -o "$sorted.bam" -O bam -

samtools index "$sorted.bam"


################################## Make Bigwig using deeptools () ###################
#using the sorted bam output, make a bigwig file
#https://deeptools.readthedocs.io/en/develop/

#create bw
bamCoverage -p $THREADS --MNase -bs 1 --smoothLength 25 -of bigwig -b "$sorted.bam" -o "$bigwig"
done
