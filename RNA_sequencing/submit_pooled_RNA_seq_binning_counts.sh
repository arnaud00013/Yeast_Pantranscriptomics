#!/bin/bash
#SBATCH --account=def-XXXXXX          # <- customize
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --array=1-6                   # one task per YPD library in list_of_samples.tsv
#SBATCH --job-name=pooled_RNAseq_binning
#
# Pre-processing, species binning (demultiplexing) and per-species gene counts
# for the six pooled YPD RNA-seq libraries (3 pools x 2 biological replicates).
# The fastp, STAR, samtools and featureCounts commands and parameters are the ones
# recorded in Yeast_subphylum_pan_transcriptomic_study_sample_preprocessing_and_QC.html.
# Usernames, module versions, paths and resources must be customized to your
# own computing environment.
#
# Prerequisite: sbatch make_multispecies_references.sh
# Usage:        sbatch submit_pooled_RNA_seq_binning_counts.sh

set -euo pipefail

module load StdEnv/2023 fastqc fastp star samtools subread   # <- set versions for your cluster

RAW_DIR=raw_RNA_seq_pools                        # FASTQ downloaded from NCBI SRA
FASTP_files_path=fastp_preprocessed
REFERENCE_FASTA_DIR=references/fasta
REFERENCE_GTF_DIR=references/gtf
mkdir -p ${FASTP_files_path} STAR_BAM read_counts_summary fastqc_raw fastqc_trimmed

M1="Ambrosiozyma_oregonensis Candida_ulmi Hanseniaspora_opuntiae Lipomyces_kononenkoae Meyerozyma_caribbica Saccharomycopsis_amapae Schwanniomyces_capriottii"
M2="Blastobotrys_peoriensis Candida_blattariae Cyberlindnera_japonica Hanseniaspora_nectarophila Nakazawaea_wyomingensis Wickerhamomyces_anomalus Wickerhamomyces_strasburgensis"
M3="Candida_dosseyi Cyberlindnera_mrakii Kodamaea_ohmeri Nakazawaea_peltata Saturnispora_besseyi Schwanniomyces_polymorphus Sugiyamaella_smithiae"

# Pick this task's sample (YPD libraries only)
line=$(tail -n +2 list_of_samples.tsv | awk -F'\t' '$10=="TRUE"' | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(echo "$line" | cut -f1)
fq=$(echo "$line" | cut -f2)
pool=$(echo "$line" | cut -f3)
SPECIES=${!pool}

#Step 1: Quality report on raw reads
fastqc -t 4 -o fastqc_raw ${RAW_DIR}/${fq}_R1_001.fastq.gz ${RAW_DIR}/${fq}_R2_001.fastq.gz

#Step 2: Clean raw reads with fastp
fastp -i ${RAW_DIR}/${fq}_R1_001.fastq.gz \
      -I ${RAW_DIR}/${fq}_R2_001.fastq.gz \
      -o ${FASTP_files_path}/${sample}_trimmed_1.fastq.gz \
      -O ${FASTP_files_path}/${sample}_trimmed_2.fastq.gz \
      -h ${FASTP_files_path}/fastp_${sample}.html \
      -j ${FASTP_files_path}/fastp_${sample}.json \
      -l 151 \
      --trim_poly_g \
      --trim_poly_a \
      --cut_front \
      --cut_front_mean_quality 20 \
      --cut_tail \
      --cut_tail_mean_quality 20 \
      --detect_adapter_for_pe \
      --overrepresentation_analysis \
      --thread 2

#Step 3: Quality report on trimmed reads
fastqc -t 4 -o fastqc_trimmed ${FASTP_files_path}/${sample}_trimmed_1.fastq.gz ${FASTP_files_path}/${sample}_trimmed_2.fastq.gz

#Step 4: Build STAR index for the sample-specific multi-species reference
STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ${sample}_star_index \
     --genomeFastaFiles ${REFERENCE_FASTA_DIR}/${sample}_multispecies.fas \
     --sjdbGTFfile ${REFERENCE_GTF_DIR}/${sample}_multispecies.gtf \
     --genomeSAindexNbases 12 \
     --sjdbOverhang 150

#Step 5: Align reads with STAR and index output bam file
STAR --genomeDir ${sample}_star_index \
     --runThreadN 4 \
     --readFilesIn ${FASTP_files_path}/${sample}_trimmed_1.fastq.gz ${FASTP_files_path}/${sample}_trimmed_2.fastq.gz \
     --readFilesCommand zcat \
     --outFilterMultimapNmax 999 \
     --outFileNamePrefix STAR_BAM/${sample}_ \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI NM MD AS \
     --outReadsUnmapped Fastx
samtools index STAR_BAM/${sample}_Aligned.sortedByCoord.out.bam

#Step 6: Split the BAM file per species and count reads per gene with featureCounts
for sp in $SPECIES; do
    # Extract contig names for this species
    grep "^>${sp}_" ${REFERENCE_FASTA_DIR}/${sample}_multispecies.fas | cut -d' ' -f1 | sed 's/>//' > ${sample}_${sp}_contigs.txt

    # Extract alignments for this species
    samtools view -b -o STAR_BAM/${sample}_${sp}.bam STAR_BAM/${sample}_Aligned.sortedByCoord.out.bam $(cat ${sample}_${sp}_contigs.txt)

    # Count reads per gene using featureCounts
    featureCounts -T 4 -p --countReadPairs -B -C \
        -a ${REFERENCE_GTF_DIR}/${sample}_${sp}_modified.gtf \
        -o read_counts_summary/${sample}_${sp}_counts.txt \
        STAR_BAM/${sample}_${sp}.bam
done

# STAR-unmapped reads (STAR_BAM/${sample}_Unmapped.out.mate1/2) are profiled
# afterwards with CCMetagen/kma against a database of the 21 species (binning QC).
# After all array tasks finish, aggregate the QC reports with: multiqc fastqc_raw fastqc_trimmed ${FASTP_files_path}
