#!/bin/bash
#This is a reduced version of the pipeline that just makes a new consensus
#Pavitra Roychoudhury, Nov 2017

#Usage: 
# 		prokka-genbank_to_fasta_db ./refs/NC_001806.2.gb ./refs/NC_001798.2.gb ./refs/KF781518.1.gb > ./refs/HSV_proteins.faa
#For paired-end library
#		hsv1_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz
#For single-end library
#		hsv1_pipeline.sh -s yourreads.fastq.gz
#This is meant to be run on the cluster (typically through sbatch) so if run locally,
#first set the environment variable manually, e.g.
#		SLURM_CPUS_PER_TASK=8
#or whatever is the number of available processors 
 
#Load required tools
# module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0-fh1
# module load prokka/1.11-foss-2016b-BioPerl-1.7.0

PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.9.0-Linux/bin:$HOME/mugsy_x86-64-v1r2.2:$HOME/last759/:$HOME/bbmap/:$HOME/samtools-1.3.1/:
export MUGSY_INSTALL=$HOME/mugsy_x86-64-v1r2.2
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK
# echo "Path: "$PATH

while getopts ":1:2:s:f" opt; do
	case $opt in
		1) in_fastq_r1="$OPTARG"
			paired="true"
		;;
		2) in_fastq_r2="$OPTARG"
			paired="true"
		;;
		s) in_fastq="$OPTARG"
			paired="false"
		;;
		f) filter="true"
		;;
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

printf "Input arguments:\n\n"
echo $@



##  PAIRED-END  ##
if [[ $paired == "true" ]]
then
if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq_r1%%_R1_001.fastq*})



##  SINGLE-END  ##
else 
if [[ $paired == "false" ]]
then
if [ -z $in_fastq ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq%%_R1_001.fastq*})


fi
fi


#Call R script to generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs_all
mkdir -p ./stats
if [[ $paired == "true" ]]
then
Rscript --vanilla hsv_generate_consensus.R s1=\"$in_fastq_r1\"
else
if [[ $paired == "false" ]]
then
Rscript --vanilla hsv_generate_consensus.R s1=\"$in_fastq\"
fi
fi

#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka_hsv1
prokka --outdir './annotations_prokka_hsv1/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 1' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv1/'$sampname/*.fa
mkdir -p ./annotations_prokka_hsv2sd90e
prokka --outdir './annotations_prokka_hsv2sd90e/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 2' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv2sd90e/'$sampname/*.fa
mkdir -p ./annotations_prokka_hsv2hg52
prokka --outdir './annotations_prokka_hsv2hg52/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 2' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv2hg52/'$sampname/*.fa
