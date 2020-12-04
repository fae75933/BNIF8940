#!/bin/bash
#SBATCH --job-name=chipseq
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/fae75933/log.%j
#SBATCH --mail-user=fae75933@uga.edu
#SBATCH --mail-type=END,FAIL

# Iterate over fastq files in a directory
# Step2 - map fastq files to N. crassa genome, output aligned reads in bam format, then sorted and indexed bam files
#Step 3 - create a bigwig file for viewing in a genome browser and for downstream analysis;

#User input needed!!! Add path to directory containing the fastq files. Include a wild card symbol at end so the script will analyze all files
THREADS=6

#Directory to iterate over with a * at the end
FILES="/scratch/fae75933/finalproject/Felicia/chipseq/*.fastq" #Don't forget the *

##manually create a directory to store output and then put the path to that output directory here for writing

OUTDIR="/scratch/fae75933/finalproject"


##make output directories that you need. These should be modified to match the software in your specific pipeline

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
        #Examples to Get Different parts of the file name
        #See here for details: http://tldp.org/LDP/abs/html/refcards.html#AEN22664
        #if you need to extract the directory but not the file, use the syntax below for variable dir
        #dir=${f%/*}    #


#create file name variables to use in the downstream analysis

#use sed to get the second read matching the input file

#filename variable for the sorted bam file
sorted="$OUTDIR/SortedBamFiles/$file"

#filename variable for the deeptools big wig output
bigwig="$OUTDIR/Bigwigs/$file.bw"


#############    Map Reads and convert to sorted bam files #########################
#http://bio-bwa.sourceforge.net/bwa.shtml
#http://www.htslib.org/doc/1.2/samtools.html

####load modules just before use, because these modules load all dependencies


##map reads and convert to sorted bam file. This is a bwa command, then output is piped to "samtools view", them this output is piped to "samtools sort"
bwa mem -M -v 3 -t $THREADS /scratch/fae75933/genomesfolder/GCA_000182925_neurospora $f | samtools view -bhu - | samtools sort -@ $THREADS -T $file -o "$sorted.bam" -O bam -

samtools index "$sorted.bam"


################################## Make Bigwig using deeptools () ###################
#using the sorted bam output, make a bigwig file
#https://deeptools.readthedocs.io/en/develop/


#create bw
bamCoverage -p $THREADS -bs 1 --smoothLength 25 -of bigwig -b "${sorted}.bam" -o "$bigwig"
#######For CutandRun data, you need to analyze the bam using the --MNase option. This
#bamCoverage -p 12 --MNase -bs 1 --smoothLength 25 -of bigwig -b "$sorted.bam" -o "$bigwig"
done
