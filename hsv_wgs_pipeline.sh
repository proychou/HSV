#!/bin/bash
#Pipeline for whole genome HSVs
#May 10, 2017

#Usage: 
#First build reference for bowtie and make a copy of the ref seq:
#		module load bowtie2
#		bowtie2-build './NC_001806.2.fasta' hsv1_ref
#		cp './NC_001806.2.fasta' hsv1_ref.fasta
#For paired-end library
#		hsv1_pipeline.sh -1 yourreads_r1.fastq.gz -2 yourreads_r2.fastq.gz -r hsv1_ref
#For single-end library
#		hsv1_pipeline.sh -s yourreads.fastq.gz -r hsv1_ref
#This is meant to be run on the cluster (typically through sbatch) so if run locally,
#first set the environment variable manually, e.g.
#		SLURM_CPUS_PER_TASK=8
#or whatever is the number of processors available for the process
 
#Load required tools
#Note that samtools, mugsy, spades, bbmap and last are all locally installed and need to be updated manually as required
# module load bowtie2
# module load FastQC/0.11.5-Java-1.8.0_92
# ml R/3.3.3-foss-2016b-fh2


PATH=$PATH:$HOME/.local/bin:$HOME/SPAdes-3.9.0-Linux/bin:$HOME/mugsy_x86-64-v1r2.2:$HOME/last759/:$HOME/bbmap/:$HOME/samtools-1.3.1/:
export MUGSY_INSTALL=$HOME/mugsy_x86-64-v1r2.2
export PATH=$PATH:$EBROOTPROKKA/bin:$EBROOTPROKKA/db:
echo "Number of cores used: "$SLURM_CPUS_PER_TASK
echo "Path: "$PATH

while getopts ":1:2:s:r:" opt; do
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
		\?) echo "Invalid option -$OPTARG" >&2
    	;;
	esac
done

#For testing single-end
# in_fastq='/fh/fast/jerome_k/HSV2_WGS/fastq_files/fastq_files_Sep2016/2003-27026-1099CUL_S300_L001_R1_001.fastq'
# paired=false


##  PAIRED-END  ##
if [[ $paired == "true" ]]
then
if [ -z $in_fastq_r1 ] || [ -z $in_fastq_r2 ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq_r1%%_R1_001.fastq*})

#FastQC report on raw reads
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq_r1 $in_fastq_r2 -o ./fastqc_reports_raw

#Adapter trimming with bbduk
mkdir -p ./trimmed_fastq
bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2  out1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'  out1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' './trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'

#Quality trimming
mkdir -p ./preprocessed_fastq
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' out1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' out2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#FastQC report on processed reads
mkdir -p ./fastqc_reports_trimmed
fastqc './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o ./fastqc_reports_trimmed

#Map reads to reference
mkdir -p ./mapped_reads
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x $ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes 
mkdir -p './contigs/'$sampname
spades.py -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

#Delete some spades folders to free up space
rm -r './contigs/'$sampname'/corrected' 



##  SINGLE-END  ##
else 
if [[ $paired == "false" ]]
then
if [ -z $in_fastq ]
then
echo "Missing input argument."
fi

sampname=$(basename ${in_fastq%%_R1_001.fastq*})

#FastQC report on raw reads
mkdir -p ./fastqc_reports_raw
fastqc $in_fastq -o ./fastqc_reports_raw

#Adapter trimming with bbduk
mkdir -p ./trimmed_fastq
bbduk.sh in=$in_fastq out='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'  out='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' ref=~/bbmap/resources/adapters.fa k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' 

#Quality trimming
mkdir -p ./preprocessed_fastq
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' out='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#FastQC report on processed reads
mkdir -p ./fastqc_reports_trimmed
fastqc './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o ./fastqc_reports_trimmed

#Map reads to reference
mkdir -p ./mapped_reads
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x $ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes
mkdir -p './contigs/'$sampname
spades.py -s './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

fi
fi



#Map reads to refs (generate sorted bam)
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
if [ -f './mapped_reads/'$sampname'_'$ref'.sam' ]
then
~/samtools-1.3.1/samtools view -bh -o './mapped_reads/'$sampname'_'$ref'.bam' './mapped_reads/'$sampname'_'$ref'.sam' -T $ref'.fasta'  
rm './mapped_reads/'$sampname'_'$ref'.sam'
~/samtools-1.3.1/samtools sort -o './mapped_reads/'$sampname'_'$ref'.sorted.bam' './mapped_reads/'$sampname'_'$ref'.bam' 
rm './mapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'Mapping to '$ref 'failed. No sam file found'
fi
done


#Map contigs to refs
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
mugsy --directory `readlink -f './contigs/'$sampname` --prefix 'aligned_scaffolds_'$ref $ref'.fasta' `readlink -f './contigs/'$sampname'/scaffolds.fasta'`
sed '/^a score=0/,$d' './contigs/'$sampname'/aligned_scaffolds_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf'
python ~/last-759/scripts/maf-convert sam -d './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
~/samtools-1.3.1/samtools view -bSq 2 -T $ref'.fasta' './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam' | ~/samtools-1.3.1/samtools sort > './contigs/'$sampname'/aligned_scaffolds_'$ref'.bam'
rm './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
done
rm *.mugsy.log


#Merge and reheader bams
mkdir -p ./merged_bam_map_plus_assembly
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
~/samtools-1.3.1/samtools view -h -o './contigs/'$sampname'/temp.sam' './mapped_reads/'$sampname'_'$ref'.sorted.bam' #for the header
~/samtools-1.3.1/samtools reheader './contigs/'$sampname'/temp.sam' './contigs/'$sampname'/aligned_scaffolds_'$ref'.bam' > './contigs/'$sampname'/temp.bam'
rm './contigs/'$sampname'/temp.sam'
~/samtools-1.3.1/samtools merge -fp './contigs/'$sampname'/merged_alignment_'$ref'.bam'  './contigs/'$sampname'/temp.bam' './mapped_reads/'$sampname'_'$ref'.sorted.bam'
rm './contigs/'$sampname'/temp.bam'
~/samtools-1.3.1/samtools sort -o './merged_bam_map_plus_assembly/'$sampname'_'$ref'.bam' './contigs/'$sampname'/merged_alignment_'$ref'.bam'
done


#Call R script to merge bams and generate a consensus sequence
mkdir -p ./consensus_seqs_all
mkdir -p ./cleaned_consensus
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
mkdir -p ./annotations_prokka
prokka --outdir './annotations_prokka/'$sampname'/' --force --kingdom 'Viruses' --genus 'Simplexvirus' --species '' --proteins HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka/'$sampname/*.fa

