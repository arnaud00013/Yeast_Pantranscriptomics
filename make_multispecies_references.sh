#!/bin/bash
#SBATCH --account=def-XXXXXX          # <- customize
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=make_multispecies_refs
#
# Build the pool-specific multi-species reference files used for species binning.
# Every contig name in each species' genome FASTA and GTF is prefixed with the
# species name ("<Genus>_<species>_<contig>") and the seven species of a pool are
# concatenated into one FASTA/GTF per sample. A prefixed single-species GTF
# (<sample>_<species>_modified.gtf) is also written for featureCounts.
#
# Reference genomes and annotations: Opulente et al. (2024) Science,
# doi:10.1126/science.adj4503 (figshare collection 6714042).
# Expected input names (customize): ${GENOME_DIR}/<Genus>_<species>.fas and
#                                   ${GENOME_DIR}/<Genus>_<species>.gtf
#
# Usage: sbatch make_multispecies_references.sh

set -euo pipefail

GENOME_DIR=/path/to/Opulente2024_genomes        # <- customize
REFERENCE_FASTA_DIR=references/fasta
REFERENCE_GTF_DIR=references/gtf
mkdir -p ${REFERENCE_FASTA_DIR} ${REFERENCE_GTF_DIR}

M1="Ambrosiozyma_oregonensis Candida_ulmi Hanseniaspora_opuntiae Lipomyces_kononenkoae Meyerozyma_caribbica Saccharomycopsis_amapae Schwanniomyces_capriottii"
M2="Blastobotrys_peoriensis Candida_blattariae Cyberlindnera_japonica Hanseniaspora_nectarophila Nakazawaea_wyomingensis Wickerhamomyces_anomalus Wickerhamomyces_strasburgensis"
M3="Candida_dosseyi Cyberlindnera_mrakii Kodamaea_ohmeri Nakazawaea_peltata Saturnispora_besseyi Schwanniomyces_polymorphus Sugiyamaella_smithiae"

for sample in $(tail -n +2 list_of_samples.tsv | awk -F'\t' '$10=="TRUE"{print $1}'); do
    pool=${sample%%_*}
    SPECIES=${!pool}
    : > ${REFERENCE_FASTA_DIR}/${sample}_multispecies.fas
    : > ${REFERENCE_GTF_DIR}/${sample}_multispecies.gtf
    for sp in $SPECIES; do
        # FASTA: ">contig ..." -> ">${sp}_contig ..."
        sed "s/^>/>${sp}_/" ${GENOME_DIR}/${sp}.fas >> ${REFERENCE_FASTA_DIR}/${sample}_multispecies.fas
        # GTF: prefix column 1 (seqname), keep comment lines untouched
        awk -F'\t' -v p="${sp}_" 'BEGIN{OFS="\t"} /^#/{print; next} {$1=p $1; print}' ${GENOME_DIR}/${sp}.gtf \
            > ${REFERENCE_GTF_DIR}/${sample}_${sp}_modified.gtf
        grep -v "^#" ${REFERENCE_GTF_DIR}/${sample}_${sp}_modified.gtf >> ${REFERENCE_GTF_DIR}/${sample}_multispecies.gtf
    done
done
