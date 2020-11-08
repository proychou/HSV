#!/bin/bash
#Pipeline for whole genome HSVs
#October 2017
#Pavitra Roychoudhury

#Usage: 
#First build reference for bowtie and make a copy of the ref seqs:
#		module load bowtie2
# 		bowtie2-build './refs/NC_001806.2.fasta' ./refs/hsv1_ref
# 		bowtie2-build './refs/NC_001798.2.fasta' ./refs/hsv2_ref_hg52
# 		bowtie2-build './refs/KF781518.1.fasta' ./refs/hsv2_sd90e
# 		cp './refs/NC_001806.2.fasta' ./refs/hsv1_ref.fasta
# 		cp './refs/NC_001798.2.fasta' ./refs/hsv2_ref_hg52.fasta
# 		cp './refs/KF781518.1.fasta' ./refs/hsv2_sd90e.fasta
# 		bwa index ./refs/hsv1_ref.fasta
# 		bwa index ./refs/hsv2_ref_hg52.fasta
# 		bwa index ./refs/hsv2_sd90e.fasta

# 		cat ./refs/hsv1_ref.fasta ./refs/hsv2_ref_hg52.fasta ./refs/hsv2_sd90e.fasta > ./refs/hsv_refs.fasta
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
#Note that samtools, mugsy, spades, bbmap and last are all locally installed and need to be updated manually as required
# module load bowtie2
# module load FastQC/0.11.5-Java-1.8.0_92
# module load R-bundle-Bioconductor/3.5-foss-2016b-R-3.4.0-fh1
# module load prokka/1.11-foss-2016b-BioPerl-1.7.0

#To do: 
# - replace local version of samtools with module load since these have now caught up
# - add a restart option
# - move some of the common parts of this script (between HHV6 and HSV) to ViralWGS 

PATH=$PATH:$HOME/.local/bin:$HOME/mugsy_x86-64-v1r2.2:$HOME/last759/:
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

#FastQC report on raw reads
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq_r1
fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq_r2 

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in1=$in_fastq_r1 in2=$in_fastq_r2  out1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz' ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'  out1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' out2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 tpe tbo overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_r1_tmp.fastq.gz' './trimmed_fastq/'$sampname'_trimmed_r2_tmp.fastq.gz'

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in1='./trimmed_fastq/'$sampname'_trimmed_r1.fastq.gz' in2='./trimmed_fastq/'$sampname'_trimmed_r2.fastq.gz' out1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' out2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#Use bbduk to filter reads that match HSV genomes
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hsv_refs.fasta ... \n\n\n"
mkdir -p ./filtered_fastq/

bbduk.sh in1='./preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' in2='./preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' out1='./filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' out2='./filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' outm1='./filtered_fastq/'$sampname'_matched_r1.fastq.gz' outm2='./filtered_fastq/'$sampname'_matched_r2.fastq.gz' ref='./refs/hsv_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hsv.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK

rm './filtered_fastq/'$sampname'_unmatched_r1.fastq.gz' './filtered_fastq/'$sampname'_unmatched_r2.fastq.gz' 
mv './filtered_fastq/'$sampname'_matched_r1.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz'
mv './filtered_fastq/'$sampname'_matched_r2.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz'
fi

#FastQC report on processed reads
mkdir -p ./fastqc_reports_preprocessed
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' 
fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' 

#Map reads to reference
printf "\n\nMapping reads to reference seqs hsv1_ref, hsv2_ref_hg52 and hsv2_sd90e ... \n\n\n"
mkdir -p ./mapped_reads
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x ./refs/$ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes 
printf "\n\nStarting de novo assembly ... \n\n\n"
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
printf "\n\nFastQC report on raw reads ... \n\n\n"
mkdir -p ./fastqc_reports_raw
fastqc -o ./fastqc_reports_raw -t $SLURM_CPUS_PER_TASK $in_fastq 

#Adapter trimming with bbduk
printf "\n\nAdapter trimming ... \n\n\n"
mkdir -p ./trimmed_fastq
bbduk.sh in=$in_fastq out='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz' ref=adapters,artifacts k=21 ktrim=r mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'  out='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' ref=adapters,artifacts k=21 ktrim=l mink=4 hdist=2 overwrite=TRUE t=$SLURM_CPUS_PER_TASK 
rm './trimmed_fastq/'$sampname'_trimmed_tmp.fastq.gz'

#Quality trimming
printf "\n\nQuality trimming ... \n\n\n"
mkdir -p ./preprocessed_fastq
bbduk.sh in='./trimmed_fastq/'$sampname'_trimmed.fastq.gz' out='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' t=$SLURM_CPUS_PER_TASK qtrim=rl trimq=20 maq=10 overwrite=TRUE minlen=20

#Use bbduk to filter reads that match HSV genomes
if [[ $filter == "true" ]]
then
printf "\n\nK-mer filtering using hsv_refs.fasta ... \n\n\n"
mkdir -p ./filtered_fastq/
bbduk.sh in='./preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' out='./filtered_fastq/'$sampname'_unmatched.fastq.gz' outm='./filtered_fastq/'$sampname'_matched.fastq.gz' ref='./refs/hsv_refs.fasta' k=31 hdist=2 stats='./filtered_fastq/'$sampname'_stats_hsv.txt' overwrite=TRUE t=$SLURM_CPUS_PER_TASK
rm './filtered_fastq/'$sampname'_unmatched.fastq.gz' 
mv './filtered_fastq/'$sampname'_matched.fastq.gz' './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz'
fi

#FastQC report on processed reads
printf "\n\nFastQC report on preprocessed reads ... \n\n\n"
mkdir -p ./fastqc_reports_preprocessed
fastqc -o ./fastqc_reports_preprocessed -t $SLURM_CPUS_PER_TASK './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' 


#Map reads to reference
printf "\n\nMapping reads to reference seqs hsv1_ref, hsv2_ref_hg52 and hsv2_sd90e ... \n\n\n"
mkdir -p ./mapped_reads
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x ./refs/$ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} -S './mapped_reads/'$sampname'_'$ref'.sam'
done

#Assemble with SPAdes
printf "\n\nStarting de novo assembly ... \n\n\n"
mkdir -p './contigs/'$sampname
spades.py -s './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -o './contigs/'$sampname --careful -t ${SLURM_CPUS_PER_TASK}

fi
fi



#Generate sorted bams for mapped reads
printf "\n\nMaking and sorting bam files ... \n\n\n"
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
if [ -f './mapped_reads/'$sampname'_'$ref'.sam' ]
then
samtools view -bh -o './mapped_reads/'$sampname'_'$ref'.bam' './mapped_reads/'$sampname'_'$ref'.sam' -T ./refs/$ref'.fasta'  
rm './mapped_reads/'$sampname'_'$ref'.sam'
samtools sort -o './mapped_reads/'$sampname'_'$ref'.sorted.bam' './mapped_reads/'$sampname'_'$ref'.bam' 
rm './mapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'Mapping to '$ref 'failed. No sam file found'
fi
done


printf "\n\nMapping scaffolds to reference seqs hsv1_ref, hsv2_ref_hg52 and hsv2_sd90e ... \n\n\n"
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
mugsy --directory `readlink -f './contigs/'$sampname` --prefix 'aligned_scaffolds_'$ref ./refs/$ref'.fasta' `readlink -f './contigs/'$sampname'/scaffolds.fasta'`
sed '/^a score=0/,$d' './contigs/'$sampname'/aligned_scaffolds_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf'
maf-convert sam -d './contigs/'$sampname'/aligned_scaffolds_nonzero_'$ref'.maf' > './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
samtools view -bS -T ./refs/$ref'.fasta' './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam' | samtools sort > './contigs/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
rm './contigs/'$sampname'/aligned_scaffolds_'$ref'.sam'
done
rm *.mugsy.log


#Make new reference sequence using scaffolds
printf "\n\nMaking a reference sequence for remapping ... \n\n\n"
mkdir -p ./ref_for_remapping
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bamfname='./contigs/'$sampname'/'$sampname'_aligned_scaffolds_'$ref'.bam'
reffname=./refs/$ref'.fasta'
Rscript --vanilla hsv_make_reference.R bamfname=\"$bamfname\" reffname=\"$reffname\" 
done


#Remap reads to "new" reference
printf "\n\nRe-mapping reads to assembled sequence ... \n\n\n"
mkdir -p ./remapped_reads

for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2-build './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref'_consensus.fasta' './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref
done


if [[ $paired == "false" ]]
then
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -U './preprocessed_fastq/'$sampname'_preprocessed.fastq.gz' -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > './remapped_reads/'$sampname'_'$ref'.bam'
done
fi
 
if [[ $paired == "true" ]]
then
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
bowtie2 -x './ref_for_remapping/'$sampname'_aligned_scaffolds_'$ref -1 './preprocessed_fastq/'$sampname'_preprocessed_paired_r1.fastq.gz' -2 './preprocessed_fastq/'$sampname'_preprocessed_paired_r2.fastq.gz' -p ${SLURM_CPUS_PER_TASK} | samtools view -bS - > './remapped_reads/'$sampname'_'$ref'.bam'
done
fi

#Make sorted bam
for ref in hsv1_ref hsv2_ref_hg52 hsv2_sd90e; do
if [ -f './remapped_reads/'$sampname'_'$ref'.bam' ]
then
samtools sort -o './remapped_reads/'$sampname'_'$ref'.sorted.bam' './remapped_reads/'$sampname'_'$ref'.bam'
rm './remapped_reads/'$sampname'_'$ref'.bam'
mv './remapped_reads/'$sampname'_'$ref'.sorted.bam'  './remapped_reads/'$sampname'_'$ref'.bam' 
else
echo 'No bam file found'
fi
done

#Call R script to merge bams and generate a consensus sequence
printf "\n\nGenerating consensus sequence ... \n\n\n"
mkdir -p ./consensus_seqs_all
mkdir -p ./stats
Rscript --vanilla hsv_generate_consensus.R sampname=\"$sampname\"

#Annotate
printf "\n\nAnnotating with prokka ... \n\n\n"
mkdir -p ./annotations_prokka_hsv1
prokka --outdir './annotations_prokka_hsv1/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 1' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv1/'$sampname/*.fa
mkdir -p ./annotations_prokka_hsv2sd90e
prokka --outdir './annotations_prokka_hsv2sd90e/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 2' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv2sd90e/'$sampname/*.fa
mkdir -p ./annotations_prokka_hsv2hg52
prokka --outdir './annotations_prokka_hsv2hg52/'$sampname'/' --force --kingdom 'Viruses' --genus 'Human herpesvirus 2' --species '' --proteins ./refs/HSV_proteins.faa --locustag '' --strain $sampname --prefix $sampname --gcode 1 --evalue 1e-9 './annotations_prokka_hsv2hg52/'$sampname/*.fa

#Clean up some files
rm './ref_for_remapping/'$sampname*'.fai'
rm './ref_for_remapping/'$sampname*'.bt2'
