#!/bin/bash
#SBATCH --job-name=finalproject
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=24gb
#SBATCH --time=24:00:00
#SBATCH --output=/work/gene8940/fae75933/log.%j
#SBATCH --mail-user=fae75933@uga.edu
#SBATCH --mail-type=END,FAIL

OUTDIR="/work/gene8940/fae75933/finalproject"

if [ ! -d $OUTDIR ]
then
    mkdir -p $OUTDIR
fi
cd $OUTDIR

Threads=6
MEMORY=24

# Use to login into sapelo
ssh fae75933@sap2test.gacrc.uga.edu

#Change into my directory
cd /scratch/fae75933

#Copy CUT&RUN data from Zack's directory into mine
cp -R /scratch/zlewis/Felicia   ./

#Create a genome folder
/c/Users/Owner/genomefolder/neurospora/GenbankAssembly12

#Download the N. crassa genome from Genbank (genome accession # GCA_000182925.2 )
curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/925/GCF_000182925.2_NC12/GCF_000182925.2_NC12_genomic.fna.gz | gunzip -c > GCA_000182925_neurospora.fna
curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.gff.gz | gunzip -c > GCA_000182925A_neurospora.fna
curl -s ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/182/925/GCA_000182925.2_NC12/GCA_000182925.2_NC12_genomic.gtf.gz | gunzip -c > GCA_000182925A2_neurospora.fna
#Use bwa index to create an indexed genome
ml BWA/0.7.17-GCC-8.3.0
bwa index GCA_000182925_neurospora.fna # if this does not work, try adding complete path
bwa index GCA_000182925A_neurospora.fna
bwa index GCA_000182925A2_neurospora.fna
#Download control ChIP-seq use the SRA toolkit "fast-dump" command
ml SRA-Toolkit/2.9.6-1-centos_linux64
#Download H3K4me2
prefetch -O /scratch/fae75933/finalproject/mapping SRR12229307
fastq-dump -I --split-files SRR12229307

#Download H3K27me/3
prefetch -O /scratch/fae75933/finalproject/mapping SRR11266616
fastq-dump -I --split-files SRR11266616

#Download input
prefetch -O /scratch/fae75933/finalproject/mapping SRR12229313
fastq-dump -I --split-files SRR12229313

#Code to move files from one spot to another
mv SRR* /scratch/fae75933/finalproject/Felicia/chipseq
#Code to view job queue
ls -l /scratch/fae75933
#cd into BNIF8940 directory to submit scrpt into cluster
sbatch
