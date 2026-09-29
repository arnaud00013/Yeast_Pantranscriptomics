#@Author=Arnaud NG
#This script summarizes the yeast subphylum pooled RNA-seq

#import libraries
library("ggplot2")
library("seqinr")
library("RColorBrewer")
library("randomcoloR")
library("FD")
library("vegan")
library("gplots")
library("lmPerm")
library("ggpubr")
library("gridExtra")
library("cluster")
library("tidyr")
library("doParallel")
library("foreach")
library("ape")
library("dplyr")
library("eulerr")
#library("plotly")

#import script arguments
#workspace with all the .bam and .tab files
output_workspace <- "C:/Users/arnau/Documents/Yeast_subphylum_pantranscriptomics/"
#set seed
set.seed(1234)
nb_permutations <- 999
#genome sizes
df_Genome_sizes <- read.csv2(file = paste0(output_workspace,"lst_Genome_sizes.txt"),sep = ",",header = T,stringsAsFactors = FALSE)
v_species_Genome_size <- df_Genome_sizes$Genome_size
names(v_species_Genome_size) <- df_Genome_sizes$smpl_sp_label
v_species_to_simplified_species_lbl  <- df_Genome_sizes$smpl_sp_label
names(v_species_to_simplified_species_lbl)  <- df_Genome_sizes$Species

#Get list of samples
lst_samples <- unname(sapply(X = (sapply(X = read.csv2(file = paste0(output_workspace,"lst_samples.txt"),sep = "\t",header = FALSE,stringsAsFactors = FALSE)[1], FUN=function(x) gsub(pattern = output_workspace,replacement = "",x,fixed = TRUE))),FUN=function(x) gsub(pattern = "_preprocessed_sorted.bam",replacement = "",x,fixed = TRUE)))
lst_samples_original <- lst_samples
#calculate total number of samples
nb_samples <- length(lst_samples)
nb_samples_original <- length(lst_samples_original)
#species and their mixture
df_species_and_mixtures <- read.csv2(file = paste0(output_workspace,"table_species_and_their_mixtures.csv"),sep = ",",header = T,stringsAsFactors = FALSE)
v_species_to_mixture <- df_species_and_mixtures$mixture
names(v_species_to_mixture) <- df_species_and_mixtures$species
v_lst_species_ID <- names(v_species_to_mixture)
v_lst_species_smpl_lbl <- unname(v_species_to_simplified_species_lbl[v_lst_species_ID])

#import mapping summary table
df_sample_mapping_stats <- read.csv2(file = paste0(output_workspace,"read_counts_summary/mapping_summary.csv"),sep = ",",header = T,stringsAsFactors = FALSE)
df_sample_mapping_stats$library_size <- df_sample_mapping_stats$Nb_mapped_reads + df_sample_mapping_stats$Nb_unmapped_reads
v_sample_library_size <- df_sample_mapping_stats$library_size
names(v_sample_library_size) <- df_sample_mapping_stats$Sample

#import sample read count at history
df_sample_read_count_history <- read.csv2(file = paste0(output_workspace,"read_counts_summary/Read_count_history.csv"),sep = ",",header = T,stringsAsFactors = FALSE)

#import the species read count in each sample and add the Unmapped reads + is_species_expected_in_sample
df_species_readCount_per_sample <- read.csv2(file = paste0(output_workspace,"read_counts_summary/species_read_count_per_sample.csv"),sep = ",",header = T,stringsAsFactors = FALSE)
df_species_readCount_per_sample$species_smpl_lbl <- v_species_to_simplified_species_lbl[df_species_readCount_per_sample$Species]
df_species_readCount_per_sample$sample_mixture <- substr(x = df_species_readCount_per_sample$Sample,start = 1,stop = 2)
df_species_readCount_per_sample$species_expected_mixture <- v_species_to_mixture[df_species_readCount_per_sample$Species]
df_species_readCount_per_sample$is_species_expected_in_sample <- df_species_readCount_per_sample$sample_mixture == df_species_readCount_per_sample$species_expected_mixture
  #add sample unmapped and reads as a "species"
for (current_sample in lst_samples){
  df_species_readCount_per_sample <- rbind(df_species_readCount_per_sample,data.frame(Sample=current_sample,Species="Unmapped",Nb_reads=subset(df_sample_mapping_stats, Sample==current_sample)$Nb_unmapped_reads,species_smpl_lbl="Unmapped",species_expected_mixture="None",sample_mixture=substr(x = current_sample,start = 1,stop = 2),is_species_expected_in_sample=FALSE))
  df_species_readCount_per_sample <- rbind(df_species_readCount_per_sample,data.frame(Sample=current_sample,Species="Ambiguous (discordant pairs,\nmulti-mapping, etc)",Nb_reads=unname(v_sample_library_size[current_sample]) - sum(subset(df_species_readCount_per_sample,Sample==current_sample)$Nb_reads),species_smpl_lbl="Ambiguous (discordant pairs,\nmulti-mapping, etc)",species_expected_mixture="None",sample_mixture=substr(x = current_sample,start = 1,stop = 2),is_species_expected_in_sample=FALSE))
}

#add the proportion of the sample library size for each species
df_species_readCount_per_sample$species_proportion_of_sample_lib_size <- df_species_readCount_per_sample$Nb_reads/v_sample_library_size[df_species_readCount_per_sample$Sample]
#add the FPKM
df_species_readCount_per_sample$fpkm <- (1E9 * as.numeric(df_species_readCount_per_sample$Nb_reads))/(as.numeric(unname(v_species_Genome_size[df_species_readCount_per_sample$species_smpl_lbl]))*as.numeric(unname(v_sample_library_size[df_species_readCount_per_sample$Sample])))
df_species_readCount_per_sample$ln1p_fpkm <- log1p(df_species_readCount_per_sample$fpkm)

#Abundance distribution expected species vs non-expected species
ggplot(data = df_species_readCount_per_sample,mapping=aes(y=ln1p_fpkm,x=is_species_expected_in_sample)) + geom_boxplot(fill="white",width=0.3,outlier.shape=NA) + geom_jitter(mapping = aes(col=is_species_expected_in_sample)) + xlab("") + ylab("Abundance (ln(FPKM+1))") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="bottom",axis.text.x = element_blank(),axis.title.x = element_blank() ) + facet_wrap(~Sample, ncol=4) + stat_compare_means(method = "wilcox",paired = F)+ labs(col="Expected species?")
#ggsave(filename = "Abundance_distribution_for_expected_vs_non_expected_species.png", path=output_workspace, width = 20, height = 30, units = "cm",dpi=1200)
#ggsave(filename = "Abundance_distribution_for_expected_vs_non_expected_species.svg", path=output_workspace, width = 20, height = 30, units = "cm",device = svg)

#create palette for species with the last 2 colors representing Unmapped reads (other species?) and ambiguous reads (discordant pairs, chimeric reads, etc)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(208)
color <- sample(col_vector, length(unique(df_species_readCount_per_sample$species_smpl_lbl)))
pie(1:length(color), col=color,labels = color)
palette_species <- sample(color,length(unique(df_species_readCount_per_sample$species_smpl_lbl)))
names(palette_species) <- unique(df_species_readCount_per_sample$species_smpl_lbl)
palette_species[length(palette_species)] <- "black"
palette_species[length(palette_species)-1] <- "red"
palette_species[unname(palette_species)=="#377EB8"]="seagreen"

#Stacked barplot for Number of paired-end fragments per species (no facets)
ggplot(data = df_species_readCount_per_sample,mapping=aes(y=Nb_reads,x=Sample,group=Sample,fill=factor(species_smpl_lbl,levels=names(palette_species)),col=as.character(is_species_expected_in_sample))) + geom_col() + xlab("Sample") + ylab("Number of paired-end fragments") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="right",axis.title.x = element_text(size=12),axis.text.x = element_text(size=9, angle = 60,hjust=1)) + labs(col="Expected species?",fill="Species") + scale_color_manual(values=c("TRUE"="black","FALSE"="red")) + scale_fill_manual(values=palette_species) 
#ggsave(filename = "Nb_reads_per_species_not_controlled_for_GenomeSize_for_expected_vs_non_expected_species.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)
#ggsave(filename = "Nb_reads_per_species_not_controlled_for_GenomeSize_for_expected_vs_non_expected_species.svg", path=output_workspace,  width = 30, height = 20, units = "cm",device = svg)

ggplot(data = df_species_readCount_per_sample,mapping=aes(y=species_proportion_of_sample_lib_size,x=Sample,group=Sample,fill=factor(species_smpl_lbl,levels=names(palette_species)),col=as.character(is_species_expected_in_sample))) + geom_col() + xlab("Sample") + ylab("Proportion of the\nsample library size") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="right",axis.title.x = element_text(size=12),axis.text.x = element_text(size=9, angle = 60,hjust=1)) + labs(col="Expected species?",fill="Species") + scale_color_manual(values=c("TRUE"="black","FALSE"="red")) + scale_fill_manual(values=palette_species) + scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))
#ggsave(filename = "Proportion_of_sample_libsize_for_expected_vs_non_expected_species.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)
#ggsave(filename = "Proportion_of_sample_libsize_for_expected_vs_non_expected_species.svg", path=output_workspace, width = 30, height = 20, units = "cm",device = svg)

#Unstacked coverage barplot with facets
df_species_readCount_per_sample$condition <- substr(x = df_species_readCount_per_sample$Sample,start = 4,stop = 12)
ggplot(data = df_species_readCount_per_sample,mapping=aes(y=Nb_reads,x=factor(species_smpl_lbl,levels=names(palette_species)),fill=as.character(is_species_expected_in_sample))) + geom_col(position="dodge") + xlab("Species") + ylab("Number of paired-end fragments") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="bottom",axis.title.x = element_text(size=12),axis.text.x = element_text(size=7, angle = 60,hjust=1))+ labs(fill="Expected species?") + scale_fill_manual(values=c("TRUE"="seagreen","FALSE"="red3")) + facet_grid(vars(condition), vars(sample_mixture))
#backup old file before running this line 
#ggsave(filename = "Unstacked_Coverage_per_sample_for_expected_vs_non_expected_species.svg", path=output_workspace, width = 30, height = 20, units = "cm",device = svg)
#write.table(x=subset(df_species_readCount_per_sample,grepl(pattern = "YPD_",x = condition,fixed = T)),file = paste0(output_workspace,"YPD_Table_species_readCount_per_sample.tsv"),sep = "\t",na = "NA",row.names = F,col.names = T)
#Number of PAIRED-END FRAGMENTS per sample at each step (paired-end fragments count history)
  #dodging bars
df_sample_read_count_history$sample_mixture <- substr(x = df_sample_read_count_history$Sample,start = 1,stop = 2)
df_sample_read_count_history$condition <- substr(x = df_sample_read_count_history$Sample,start = 4,stop = 12)
ggplot(data = df_sample_read_count_history,mapping=aes(y=Count,x=factor(Dataset,levels=c("Raw data","After pre-processing","After binning")))) + geom_col(position = "dodge") + xlab("Dataset") + ylab("Number of paired-end fragments") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="right",axis.title.x = element_text(size=12),axis.text.x = element_text(size=9, angle = 60,hjust=1)) + facet_grid(vars(condition), vars(sample_mixture))#+ scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))
#ggsave(filename = "Paired_End_frags_count_history_Unstacked.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)
#ggsave(filename = "Paired_End_frags_count_history_Unstacked.svg", path=output_workspace, width = 30, height = 20, units = "cm",device = svg)

  #overlapping bars NOT STACKED
ggplot(data = df_sample_read_count_history,mapping=aes(y=Count,x=Sample,fill= factor(Dataset,levels=c("Raw data","After pre-processing","After binning")))) + geom_col(position = "identity") + xlab("Sample") + ylab("Number of paired-end fragments") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),legend.title = element_text(size=12),legend.text = element_text(size=10),legend.position="right",axis.title.x = element_text(size=12),axis.text.x = element_text(size=9, angle = 60,hjust=1)) + labs(fill="Dataset") + scale_fill_brewer(palette = "Set1") #+ scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))
#ggsave(filename = "Paired_End_frags_count_history_Overlapping_bars.png", path=output_workspace, width = 30, height = 20, units = "cm",dpi=1200)
#ggsave(filename = "Paired_End_frags_count_history_Overlapping_bars.svg", path=output_workspace, width = 30, height = 20, units = "cm",device = svg)

#Only species in samples in which they have at least 0.5M reads (avg nb reads per gene = 500) NO MATTER IF THE SPECIES ARE EXPECTED OR NOT
df_wellcovered_species_readCount_per_sample <- subset(df_species_readCount_per_sample,(!is.na(fpkm))&(Nb_reads>5E5))
df_wellcovered_species_readCount_per_sample$replicate <- substr(df_wellcovered_species_readCount_per_sample$condition,start = nchar(df_wellcovered_species_readCount_per_sample$condition),stop = nchar(df_wellcovered_species_readCount_per_sample$condition))
df_wellcovered_species_readCount_per_sample$condition_X_replicate <- df_wellcovered_species_readCount_per_sample$condition
df_wellcovered_species_readCount_per_sample$condition <- substr(x = df_wellcovered_species_readCount_per_sample$condition,start = 1,stop = 3)
df_wellcovered_species_readCount_per_sample$species_sample_filename <- paste0(df_wellcovered_species_readCount_per_sample$Sample,"_",df_wellcovered_species_readCount_per_sample$Species,"_counts.txt")
#Per species stats to evaluate the species recovery from pooled RNA-seq
  #nb samples well covered; nb conditions well covered; avg nb replicates per condition
df_per_species_recovery_stats <- df_wellcovered_species_readCount_per_sample %>%
  group_by(species_smpl_lbl) %>%
  dplyr::summarise(nb_samples_well_covered = length(unique(Sample)), nb_conditions_well_covered = length(unique(condition)) )
df_per_species_recovery_stats$avg_nb_replicates <- vapply(X = df_per_species_recovery_stats$species_smpl_lbl,FUN = function(the_current_species) length(unique(subset(df_wellcovered_species_readCount_per_sample,species_smpl_lbl==the_current_species)$condition_X_replicate))/length(unique(subset(df_wellcovered_species_readCount_per_sample,species_smpl_lbl==the_current_species)$condition)),FUN.VALUE = 0.0)
df_per_species_recovery_stats$nb_conditions_with_all_replicates <- vapply(X = df_per_species_recovery_stats$species_smpl_lbl,FUN = function(the_current_species) round(length(unique(subset(df_wellcovered_species_readCount_per_sample,species_smpl_lbl==the_current_species)$condition_X_replicate))/2),FUN.VALUE = 0)

#Per condition stats to evaluate the species recovery from pooled RNA-seq
  #nb species well covered; nb_replicates recovered
df_per_condition_recovery_stats <- df_wellcovered_species_readCount_per_sample %>%
  group_by(condition) %>%
  dplyr::summarise(nb_species_well_covered = length(unique(species_smpl_lbl)), nb_replicates_covered = length(unique(replicate)) )
df_per_condition_recovery_stats$nb_species_in_both_replicates <- vapply(X = df_per_condition_recovery_stats$condition,FUN = function(the_cond) length(intersect(subset(df_wellcovered_species_readCount_per_sample,condition_X_replicate==paste0(the_cond,"_RNA_1"))$species_smpl_lbl,subset(df_wellcovered_species_readCount_per_sample,condition_X_replicate==paste0(the_cond,"_RNA_2"))$species_smpl_lbl)),FUN.VALUE = 0)

#List of samples to analyze for the non-plastic expression part (YPD and YNB)
v_lst_samples_to_analyze_non_plastic_expression <- with(df_wellcovered_species_readCount_per_sample,paste0(Sample,"_",Species,"_counts.txt"))  #subset(,condition=="YPD")

#*******************************Prepare RNA-seq data for analysis (make sure no batch effects, create relevant matrices, match gene IDs, etc)**************************************************************************#
#species gene content and paralogs
  #raw orthogroups presence/absence (Orthofinder)
df_orthogroups_per_species <- read.csv(file = paste0(output_workspace,"22781714/y1000p_orthofinder/Orthogroups/Orthogroups.tsv"),header = T,sep = "\t",na.strings = "")
rownames(df_orthogroups_per_species) <- df_orthogroups_per_species$Orthogroup
  #create gene paralogs/orthogroup/length dataframe
v_simplified_species_lbl_to_species <- names(v_species_to_simplified_species_lbl)
names(v_simplified_species_lbl_to_species) <- v_species_to_simplified_species_lbl
v_species_smpl_lbl_to_species_colnames_in_orthofile <- vapply(X = v_lst_species_smpl_lbl,FUN = function(the_smpl_sp_lbl) colnames(df_orthogroups_per_species)[grepl(pattern = tolower(v_simplified_species_lbl_to_species[the_smpl_sp_lbl]),x = tolower(colnames(df_orthogroups_per_species)),fixed = T)],FUN.VALUE = "")
names(v_species_smpl_lbl_to_species_colnames_in_orthofile) <- v_lst_species_smpl_lbl
v_species_colnames_in_orthofile_to_species_smpl_lbl <- names(v_species_smpl_lbl_to_species_colnames_in_orthofile)
names(v_species_colnames_in_orthofile_to_species_smpl_lbl) <- v_species_smpl_lbl_to_species_colnames_in_orthofile
  #only focus on the species under study
df_orthogroups_per_species <- df_orthogroups_per_species[,colnames(df_orthogroups_per_species)%in%unname(v_species_smpl_lbl_to_species_colnames_in_orthofile)]
  #parse orthofinder results and gtf files to extract gene features
#current_iteration <- 1
#nb_iterations_max <- prod(nrow(df_orthogroups_per_species),length(v_species_smpl_lbl_to_species_colnames_in_orthofile),30)
nb_cores <- detectCores() - 6 #Use 10 cores (16-6)
cl <- makeCluster(nb_cores)
registerDoParallel(cl)

    #function to collect the species orthogroups gene metadata
get_species_genes_metadata <- function(v_lst_species_current_core){
  the_it <- 1
  nb_iters <- prod(nrow(df_orthogroups_per_species),length(v_lst_species_current_core),30)
  current_species_df_per_species_orthgroups_gene_features <- NULL
  for (current_species_colname in unname(v_lst_species_current_core)){
    current_species <- unname(v_species_colnames_in_orthofile_to_species_smpl_lbl[current_species_colname])
    #import species gtf
    df_current_species_gtf <- read.csv(file = paste0(output_workspace,"gtfs/",current_species_colname,".gtf"),header = F,sep = "\t",na.strings = "")
    colnames(df_current_species_gtf) <- c("scaffold", "source", "feature","start","end","score","strand","frame","attribute")
    for (current_orthogroup in rownames(df_orthogroups_per_species)){
      #save each gene features FOR THE CURRENT SPECIES
      if (!is.na(df_orthogroups_per_species[current_orthogroup,current_species_colname])){
        current_lst_genes <- strsplit(x = df_orthogroups_per_species[current_orthogroup,current_species_colname],split = ", ",fixed = T)[[1]]
        for (current_gene in current_lst_genes){
          current_gene_length_in_current_species <- sum(abs(subset(df_current_species_gtf,grepl(pattern = current_gene,x = attribute,fixed = T)&(feature=="CDS"))$end-subset(df_current_species_gtf,grepl(pattern = current_gene,x = attribute,fixed = T)&(feature=="CDS"))$start)+1)
          current_species_df_per_species_orthgroups_gene_features <- rbind(current_species_df_per_species_orthgroups_gene_features,data.frame(gene=current_gene,species_smpl_lbl=current_species,orthogroup=current_orthogroup,paralogs_in_species=ifelse(test = length(current_lst_genes)>1,yes =  paste0(current_lst_genes[current_lst_genes!=current_gene],collapse=","),no = NA),nb_paralogs_excluding_current_gene=ifelse(test = length(current_lst_genes)>1,yes =  length(current_lst_genes)-1,no = 0),gene_length_in_species=current_gene_length_in_current_species,stringsAsFactors = F))
          print(paste0("Iteration ",the_it," out of ",nb_iters," done for species ",current_species," (if the number of genes per orthogroup is 30 on avg)!"))
          the_it <- the_it + 1
        }
      }else{
        print(paste0("Iteration ",the_it," out of ",nb_iters," done for species ",current_species," (if the number of genes per orthogroup is 30 on avg)!"))
        the_it <- the_it + 1
        next
      }
      
    }
    print(paste0("Orthogroup ",current_orthogroup,"analysis done!"))
    print(current_orthogroup)
  }
  rownames(current_species_df_per_species_orthgroups_gene_features) <- paste0(current_species_df_per_species_orthgroups_gene_features$gene," X ",current_species_df_per_species_orthgroups_gene_features$species_smpl_lbl)
  return(current_species_df_per_species_orthgroups_gene_features)
}
  #split species in a list
list_chunks_species <- split(unname(v_species_smpl_lbl_to_species_colnames_in_orthofile), cut(seq_along(1:length(unname(v_species_smpl_lbl_to_species_colnames_in_orthofile))), nb_cores, labels = FALSE))
df_per_species_orthgroups_gene_features <- foreach(chunk = list_chunks_species, .combine = rbind) %dopar% {
  get_species_genes_metadata(chunk)
}
  #stop multicore cluster 
stopCluster(cl)

#Save df_per_species_orthgroups_gene_features object in a file
saveRDS(object = df_per_species_orthgroups_gene_features,file = paste0(output_workspace,"df_per_species_orthgroups_gene_features.rds"))

#import their read counts per gene
  #Create a mega dataframe with all the samples count using rbind, naming columns appropriately and adding columns for species, sample names, condition, replicate, replicate_x_cond, etc
    #add first species sample
current_sample <- v_lst_samples_to_analyze_non_plastic_expression[1]
df_yeasts_pantranscriptome <- read.csv(file = paste0(output_workspace,"read_counts_summary/",current_sample),header = T,sep = "\t",skip = 1)
colnames(df_yeasts_pantranscriptome)[length(colnames(df_yeasts_pantranscriptome))] <- "fragment_pairs_count"
df_yeasts_pantranscriptome$species_sample_name <- current_sample
df_yeasts_pantranscriptome$species_ID <- strsplit(current_sample,split = "_",fixed = T)[[1]][7]
df_yeasts_pantranscriptome$species_smpl_lbl <- paste(strsplit(current_sample,split = "_",fixed = T)[[1]][8],strsplit(current_sample,split = "_",fixed = T)[[1]][9])
df_yeasts_pantranscriptome$condition <- strsplit(current_sample,split = "_",fixed = T)[[1]][2]
df_yeasts_pantranscriptome$replicate <- paste0("RNA",strsplit(current_sample,split = "_",fixed = T)[[1]][4])
df_yeasts_pantranscriptome$condition_X_replicate <- paste0(df_yeasts_pantranscriptome$condition,"_",df_yeasts_pantranscriptome$replicate)
df_yeasts_pantranscriptome$mixture <- strsplit(current_sample,split = "_",fixed = T)[[1]][1]
    #add other species samples iteratively
for (current_sample in v_lst_samples_to_analyze_non_plastic_expression[2:length(v_lst_samples_to_analyze_non_plastic_expression)]){
  df_to_add <- read.csv(file = paste0(output_workspace,"read_counts_summary/",current_sample),header = T,sep = "\t",skip = 1)
  colnames(df_to_add)[length(colnames(df_to_add))] <- "fragment_pairs_count"
  df_to_add$species_sample_name <- current_sample
  df_to_add$species_ID <- strsplit(current_sample,split = "_",fixed = T)[[1]][7]
  df_to_add$species_smpl_lbl <- paste(strsplit(current_sample,split = "_",fixed = T)[[1]][8],strsplit(current_sample,split = "_",fixed = T)[[1]][9])
  df_to_add$condition <- strsplit(current_sample,split = "_",fixed = T)[[1]][2]
  df_to_add$replicate <- paste0("RNA",strsplit(current_sample,split = "_",fixed = T)[[1]][4])
  df_to_add$condition_X_replicate <- paste0(df_to_add$condition,"_",df_to_add$replicate)
  df_to_add$mixture <- strsplit(current_sample,split = "_",fixed = T)[[1]][1]
  df_yeasts_pantranscriptome <- rbind(df_yeasts_pantranscriptome, df_to_add)
}

#create species gtfs list
lst_species_gtfs <- list()
for (the_species_filename in names(v_species_colnames_in_orthofile_to_species_smpl_lbl)){
  lst_species_gtfs[[the_species_filename]] <- read.csv(file =paste0(output_workspace,"gtfs/",the_species_filename,".gtf"),header = F,sep = "\t",na.strings = "")
  colnames(lst_species_gtfs[[the_species_filename]]) <- c("scaffold", "source", "feature","start","end","score","strand","frame","attribute")
  
}
  #Add FPKM and genomic presence/absence to the pantranscriptome dataframe
colnames(df_yeasts_pantranscriptome)[colnames(df_yeasts_pantranscriptome)=="species_smpl_lbl"] <- "species_all_lowercase"
df_yeasts_pantranscriptome$species_smpl_lbl <- NA
df_yeasts_pantranscriptome$FPKM <- NA
df_yeasts_pantranscriptome$int_is_detected_in_genome <- NA
for (i in 1:nrow(df_yeasts_pantranscriptome)){
  df_yeasts_pantranscriptome$species_smpl_lbl[i] <-  paste0(toupper(substr(df_yeasts_pantranscriptome$species_all_lowercase[i],start = 1,stop = 1)),substr(df_yeasts_pantranscriptome$species_all_lowercase[i],start = 2,stop = nchar(df_yeasts_pantranscriptome$species_all_lowercase[i])) )
  current_species_filename <- unname(v_species_smpl_lbl_to_species_colnames_in_orthofile[df_yeasts_pantranscriptome$species_smpl_lbl[i]])
  df_yeasts_pantranscriptome$FPKM[i] <- as.numeric((1E9)*df_yeasts_pantranscriptome$fragment_pairs_count[i])/(as.numeric(subset(df_wellcovered_species_readCount_per_sample,species_sample_filename == df_yeasts_pantranscriptome$species_sample_name[i])$Nb_reads) * as.numeric(df_yeasts_pantranscriptome$Length[i]))
  df_yeasts_pantranscriptome$int_is_detected_in_genome[i] <- as.integer(nrow(subset(lst_species_gtfs[[current_species_filename]],grepl(pattern = df_yeasts_pantranscriptome$Geneid[i],x = lst_species_gtfs[[current_species_filename]]$attribute,fixed = T)&(feature=="CDS")))>=1)
  if (i%%100==0){
    print(i)
  }
}
df_yeasts_pantranscriptome$bool_is_detected_in_genome <- as.logical(df_yeasts_pantranscriptome$int_is_detected_in_genome)
df_yeasts_pantranscriptome$bool_is_expressed <- df_yeasts_pantranscriptome$fragment_pairs_count>0
df_yeasts_pantranscriptome$Orthogroup <- df_per_species_orthgroups_gene_features[paste0(df_yeasts_pantranscriptome$Geneid,".m1 X ",df_yeasts_pantranscriptome$species_smpl_lbl),"orthogroup"]

#remove genes that are not assigned to Orthogroups (~1%)
whole_df_yeasts_pantranscriptome <- df_yeasts_pantranscriptome
df_yeasts_pantranscriptome <- subset(df_yeasts_pantranscriptome,!is.na(Orthogroup))

#make sure all expressed genes are detected at the genomic level
print(all(df_yeasts_pantranscriptome$int_is_detected_in_genome==1))

#Get FPKM in log scale 1p
df_yeasts_pantranscriptome$log2_1p_FPKM  <- log2(df_yeasts_pantranscriptome$FPKM+1)
#Save df_yeasts_pantranscriptome object in a file
##saveRDS(object = df_yeasts_pantranscriptome,file = paste0(output_workspace,"df_yeasts_pantranscriptome.rds"))

#Create a dataframe of species pairwise phylogenetic distance and distances to root
  #pairwise phylogenetic distance for species well recovered in each condition
yeasts_fulltree <- read.tree(paste0(output_workspace,"ika1qPjBmYTXt46RHctu9A_newick.txt"))
yeasts_fulltree_dist_matrix <- cophenetic.phylo(yeasts_fulltree)
mtx_phylo_dist_the_21_sps <- yeasts_fulltree_dist_matrix[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl,fixed = T), gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl,fixed = T)]
v_lst_species_smpl_lbl_YPD <- sort(unique(subset(df_wellcovered_species_readCount_per_sample,condition=="YPD")$species_smpl_lbl ))
mtx_phylo_dist_the_21_sps_YPD <- mtx_phylo_dist_the_21_sps[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD,fixed = T), gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD,fixed = T)]
v_lst_species_smpl_lbl_YNB <- sort(unique(subset(df_wellcovered_species_readCount_per_sample,condition=="YNB")$species_smpl_lbl ))
mtx_phylo_dist_the_21_sps_YNB <- mtx_phylo_dist_the_21_sps[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YNB,fixed = T), gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YNB,fixed = T)]

  #Distance to root
root_to_tip <- node.depth.edgelength(yeasts_fulltree)
  #Extract only the tips (species)
tip_distances <- root_to_tip[1:length(yeasts_fulltree$tip.label)]
names(tip_distances) <- yeasts_fulltree$tip.label

  #Generate all combinations of size 2 and add features of interest
    #YPD
the_21_sps_combinations_list_YPD <- combn(rownames(mtx_phylo_dist_the_21_sps_YPD), 2, simplify = FALSE)
df_the_21_sps_combinations_YPD <- do.call(rbind, lapply(the_21_sps_combinations_list_YPD, function(x) data.frame(Col1 = x[1], Col2 = x[2])))
colnames(df_the_21_sps_combinations_YPD) <- c("Species1","Species2")
df_the_21_sps_combinations_YPD$species_smpl_lbl_1 <- gsub(pattern = "_",replacement = " ",x = df_the_21_sps_combinations_YPD$Species1,fixed = T)
df_the_21_sps_combinations_YPD$species_smpl_lbl_2 <- gsub(pattern = "_",replacement = " ",x = df_the_21_sps_combinations_YPD$Species2,fixed = T)
df_the_21_sps_combinations_YPD$Pairwise_phylogenetic_distance <- vapply(X = 1:nrow(df_the_21_sps_combinations_YPD),FUN = function(the_i) mtx_phylo_dist_the_21_sps_YPD[df_the_21_sps_combinations_YPD$Species1[the_i],df_the_21_sps_combinations_YPD$Species2[the_i]],FUN.VALUE = 0.0)
df_the_21_sps_combinations_YPD$distance_to_root_sp_1 <- unname(tip_distances[df_the_21_sps_combinations_YPD$Species1])
df_the_21_sps_combinations_YPD$distance_to_root_sp_2 <- unname(tip_distances[df_the_21_sps_combinations_YPD$Species2])
df_the_21_sps_combinations_YPD$int_which_species_is_closer_to_root <- vapply(X = 1:nrow(df_the_21_sps_combinations_YPD),FUN = function(the_i) which.min(x = c(df_the_21_sps_combinations_YPD$distance_to_root_sp_1[the_i],df_the_21_sps_combinations_YPD$distance_to_root_sp_2[the_i])),FUN.VALUE = 1L)
df_the_21_sps_combinations_YPD$int_which_species_emerged_most_recently <- vapply(X = 1:nrow(df_the_21_sps_combinations_YPD),FUN = function(the_i) which.max(x = c(df_the_21_sps_combinations_YPD$distance_to_root_sp_1[the_i],df_the_21_sps_combinations_YPD$distance_to_root_sp_2[the_i])),FUN.VALUE = 1L)
    #YNB
the_21_sps_combinations_list_YNB <- combn(rownames(mtx_phylo_dist_the_21_sps_YNB), 2, simplify = FALSE)
df_the_21_sps_combinations_YNB <- do.call(rbind, lapply(the_21_sps_combinations_list_YNB, function(x) data.frame(Col1 = x[1], Col2 = x[2])))
colnames(df_the_21_sps_combinations_YNB) <- c("Species1","Species2")
df_the_21_sps_combinations_YNB$species_smpl_lbl_1 <- gsub(pattern = "_",replacement = " ",x = df_the_21_sps_combinations_YNB$Species1,fixed = T)
df_the_21_sps_combinations_YNB$species_smpl_lbl_2 <- gsub(pattern = "_",replacement = " ",x = df_the_21_sps_combinations_YNB$Species2,fixed = T)
df_the_21_sps_combinations_YNB$Pairwise_phylogenetic_distance <- vapply(X = 1:nrow(df_the_21_sps_combinations_YNB),FUN = function(the_i) mtx_phylo_dist_the_21_sps_YNB[df_the_21_sps_combinations_YNB$Species1[the_i],df_the_21_sps_combinations_YNB$Species2[the_i]],FUN.VALUE = 0.0)
df_the_21_sps_combinations_YNB$distance_to_root_sp_1 <- unname(tip_distances[df_the_21_sps_combinations_YNB$Species1])
df_the_21_sps_combinations_YNB$distance_to_root_sp_2 <- unname(tip_distances[df_the_21_sps_combinations_YNB$Species2])
df_the_21_sps_combinations_YNB$int_which_species_is_closer_to_root <- vapply(X = 1:nrow(df_the_21_sps_combinations_YNB),FUN = function(the_i) which.min(x = c(df_the_21_sps_combinations_YNB$distance_to_root_sp_1[the_i],df_the_21_sps_combinations_YNB$distance_to_root_sp_2[the_i])),FUN.VALUE = 1L)
df_the_21_sps_combinations_YNB$int_which_species_emerged_most_recently <- vapply(X = 1:nrow(df_the_21_sps_combinations_YNB),FUN = function(the_i) which.max(x = c(df_the_21_sps_combinations_YNB$distance_to_root_sp_1[the_i],df_the_21_sps_combinations_YNB$distance_to_root_sp_2[the_i])),FUN.VALUE = 1L)

#import species KEGG annotations in a list
df_yeasts_pantranscriptome$lbl_species_for_kegg_file <- paste0(gsub(pattern = " ",replacement = "_",x = df_yeasts_pantranscriptome$species_smpl_lbl,fixed = T),"_",df_yeasts_pantranscriptome$species_ID,".txt")
lst_species_KEGG_annotations <- list()
for (the_species_kegg_file in sort(unique(df_yeasts_pantranscriptome$lbl_species_for_kegg_file))){
  lst_species_KEGG_annotations[[the_species_kegg_file]] <- read.csv(file =paste0(output_workspace,"lst_species_genes_Kegg_Orthogroups/",the_species_kegg_file),header = F,sep = "\t",na.strings = "")
  colnames(lst_species_KEGG_annotations[[the_species_kegg_file]]) <- c("Geneid","KEGG_annotation")
}
  #Add KEGG Orthology IDs to df_yeasts_pantranscriptome
df_yeasts_pantranscriptome$Kegg_Orthorlogy_ID <- vapply(X = 1:nrow(df_yeasts_pantranscriptome),FUN = function(the_i) ifelse(test = paste0(df_yeasts_pantranscriptome$Geneid[the_i],".m1")%in%(lst_species_Kegg_Orthorlogy_IDs[[df_yeasts_pantranscriptome$lbl_species_for_kegg_file[the_i]]])$Geneid,
       yes = subset(lst_species_Kegg_Orthorlogy_IDs[[df_yeasts_pantranscriptome$lbl_species_for_kegg_file[the_i]]],Geneid==paste0(df_yeasts_pantranscriptome$Geneid[the_i],".m1") )$Kegg_Orthorlogy_ID,
       no = "NA"),FUN.VALUE = "")
df_yeasts_pantranscriptome$Kegg_Orthorlogy_ID[df_yeasts_pantranscriptome$Kegg_Orthorlogy_ID=="NA"] <- NA

#Function to get initials
get_initials <- function(name_string) {
  #Split the name string into individual words
  words <- strsplit(name_string, " ")[[1]]
  
  #Extract the first letter of each word and convert to uppercase
  initials <- (substr(words, 1, 1)) #toupper
  
  #Concatenate the initials
  paste0(paste(initials, collapse = "."),".")
}

#Cast df_yeasts_pantranscriptome to log2_FPKM matrices based on condition (best aggregate function is mean)
#sample ~ Orthogroup log2_1p_FPKM matrix
df_yeasts_pantranscriptome$short_lbl_sample <- substr(x = df_yeasts_pantranscriptome$species_sample_name,start = 1,stop = 22)
v_species_sample_name_to_species_initials <- vapply(X = sort(unique(df_yeasts_pantranscriptome$species_sample_name)),FUN = function(the_lbl) get_initials(subset(df_yeasts_pantranscriptome,species_sample_name==the_lbl)$species_smpl_lbl[1]),FUN.VALUE = "")
mtx_log2_1p_sp_sample_FPKM <- reshape2::acast(df_yeasts_pantranscriptome, species_sample_name ~ Orthogroup, value.var = "log2_1p_FPKM",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)

#save the list of Kegg Orthology IDs
v_lst_KO_IDs <- sort(unique(df_yeasts_pantranscriptome$Kegg_Orthorlogy_ID))
write.table(x=v_lst_KO_IDs,file = paste0(output_workspace,"List_expressed_genes_KO_IDs.tsv"),sep = "\t",na = "NA",row.names = F,col.names = F,quote = F)
#add KEGG pathways to df_yeasts_pantranscriptome
df_KO_ID_to_KEGG_PATHWAY <- read.delim(paste0(output_workspace,"/Table_Kegg_pathway_metadata.tsv"), header=T,sep="\t")
names(df_KO_ID_to_KEGG_PATHWAY)[1] <- "Kegg_Orthorlogy_ID"
df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY <- dplyr::left_join(x = df_yeasts_pantranscriptome,y = df_KO_ID_to_KEGG_PATHWAY,by = "Kegg_Orthorlogy_ID")
write.table(x=sort(unique(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY$Kegg_pathway)),file = paste0(output_workspace,"List_expressed_genes_KEGG_pathways.tsv"),sep = "\t",na = "NA",row.names = F,col.names = F,quote = F)
#add KEGG pathway title and KEGG class

#save df_yeasts_pantranscriptome
write.table(x=df_yeasts_pantranscriptome,file = paste0(output_workspace,"Table_df_yeasts_pantranscriptome.tsv"),sep = "\t",na = "NA",row.names = F,col.names = T)
#save df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY
write.table(x=df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,file = paste0(output_workspace,"Table_df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY.tsv"),sep = "\t",na = "NA",row.names = F,col.names = T)


#Find the best aggregate function for each orthogroup based on PCA controls:
#Verify positive controls (differences between conditions), negative controls (similarity between replicates) and batch effects (cluster by mixture?) in dimensionality reduction

  #Create a PCA of all sample FPKM (label is species initials, shape is mixture, fill color is condition and stroke color is replicate)
#Perform the PCA
sample_pca <- prcomp(mtx_log2_1p_sp_sample_FPKM)
#transform pca matrix into a dataframe and add species sample features of interest (species
#initials, mixtures, condition and replicate)
df_sample_pca <- as.data.frame(as_tibble(sample_pca$x, rownames = "species_sample"))
df_sample_pca$species_initials <- v_species_sample_name_to_species_initials[df_sample_pca$species_sample]
v_species_sample_name_to_species_smpl_lbl <- vapply(X = sort(unique(df_yeasts_pantranscriptome$species_sample_name)),FUN = function(the_lbl) subset(df_yeasts_pantranscriptome,species_sample_name==the_lbl)$species_smpl_lbl[1],FUN.VALUE = "")
df_sample_pca$species_smpl_lbl <- v_species_sample_name_to_species_smpl_lbl[df_sample_pca$species_sample]
df_sample_pca$mixture <- substr(x = df_sample_pca$species_sample,start = 1,stop = 2)
df_sample_pca$condition <- substr(x = df_sample_pca$species_sample,start = 4,stop = 6)
df_sample_pca$replicate <- substr(x = df_sample_pca$species_sample,start = 8,stop = 12)
df_sample_pca$species_initials_X_Mixture <- paste0(df_sample_pca$species_initials,"_",df_sample_pca$mixture)
  #Eigenvalues
sample_pc_eigenvalues <- sample_pca$sdev^2
df_sample_pc_eigenvalues <- tibble(PC = factor(1:length(sample_pc_eigenvalues)), 
                         variance = sample_pc_eigenvalues) %>% 
  #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
  #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

nb_pcs_recapitulating_species_sample_pca <- as.integer(subset(df_sample_pc_eigenvalues,explained_variance_cum>=0.99)$PC)[1]
  #Visualize sample PCA eigen values
df_sample_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = explained_variance)) +
  geom_line(aes(y = explained_variance_cum, group = 1)) + 
  geom_point(aes(y = explained_variance_cum)) +
  geom_hline(yintercept = 0.99,col="red",lty=2)+
  geom_vline(xintercept = nb_pcs_recapitulating_species_sample_pca,col="red",lty=2)+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  labs(x = "Principal component", y = "Fraction variance explained")
ggsave(filename = "PCA_species_sample_Cumul_Explained_variance.png", path=output_workspace, width = 25/2.54, height = 15/2.54)

  #Visualize sample PCA
df_sample_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(aes(shape = condition, col = replicate), size = 2, stroke = 2) + 
  ggrepel::geom_text_repel(aes(label = species_initials_X_Mixture)) #ggrepel::geom_text_repel
ggsave(filename = "PCA_species_sample.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

subset(df_sample_pca,condition=="YPD") %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = paste0("PC1 (",signif(df_sample_pc_eigenvalues$explained_variance[1]*100,4),"%)"), y = paste0("PC2 (",signif(df_sample_pc_eigenvalues$explained_variance[2]*100,4),"%)"))+
  geom_point(aes(col = replicate), size = 2, stroke = 2)  +
  ggrepel::geom_text_repel(aes(label = species_initials_X_Mixture)) #ggrepel::geom_text_repel
ggsave(filename = "YPD_PCA_species_sample.svg", path=output_workspace, width = 20/2.54, height = 15/2.54)
ggsave(filename = "YPD_PCA_species_sample.pdf", path=output_workspace, width = 20/2.54, height = 15/2.54,dpi=600)
ggsave(filename = "YPD_PCA_species_sample.png", path=output_workspace, width = 20/2.54, height = 15/2.54,dpi=600)

  #Visualize top 10 ORTHOGROUPS (or "gene") with highest loadings in sample PCA's pc1 and pc2
pc_loadings <- sample_pca$rotation %>% 
  as_tibble(rownames = "gene") %>%
  as.data.frame()

v_sample_pc_eigenvalues <- sample_pc_eigenvalues
names(v_sample_pc_eigenvalues) <- paste0("PC",1:length(sample_pc_eigenvalues))

#as pc1/pc2 loading dataframe
df_sorted_pc1_and_pc2_loadings <- pc_loadings %>% 
  #select only the PCs we are interested in
  dplyr::select(gene, PC1, PC2) %>%
  #convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  #for each PC
  dplyr::group_by(PC) %>%
  as.data.frame()

df_sorted_pc1_and_pc2_loadings$the_princ_comp_explained_var <- v_sample_pc_eigenvalues[df_sorted_pc1_and_pc2_loadings$PC]
df_sorted_pc1_and_pc2_loadings$normalized_loading <- df_sorted_pc1_and_pc2_loadings$loading*df_sorted_pc1_and_pc2_loadings$the_princ_comp_explained_var


df_sorted_pc1_and_pc2_loadings <- df_sorted_pc1_and_pc2_loadings %>% 
  #arrange by descending order of loading
  dplyr::arrange(desc(abs(normalized_loading))) 

top10_gene_orthogroup <- df_sorted_pc1_and_pc2_loadings$gene[1:10]

#print top 10 genes with highest loadings in sample PCA's pc1 and pc2
  #as vector
top10_gene_orthogroup
top10_gene_orthogroup_functions <- c("Cell wall organization mannoprotein CWP2", "unknown","Cell membrane channel/ Gap junction alpha-8 protein GJA8", "ribosomal protein","unknown", "unknown", "unknown", "unknown", "unknown", "unknown")
top10_gene_orthogroup_functions

#visualize sample PC1/PC2 top 10 orthogroup pc loading
top_loadings <- pc_loadings %>% 
  dplyr::filter(gene %in% top10_gene_orthogroup)

ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               col = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
            nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))+
  labs(x = "PC1", y = "PC2")
ggsave(filename = "Top10_pc1_pc2_loadings_PCA_species_sample.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

  #Quantify PCA metrics (get sample distance and compare distances distribution in function of is_mixture_the_same and is_replicate_of_same_sample OR correlate with species phylo distance)
mtx_species_sample_distance_in_PCA <- as.matrix(vegan::vegdist(sample_pca$x, method = "euclidean"))
#save the PCA stats as tsv files
  #save pca matrix
write.table(x=sample_pca$x,file = paste0(output_workspace,"mtx_species_sample_PCA.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
  #save pca disance matrix
write.table(x=mtx_species_sample_distance_in_PCA,file = paste0(output_workspace,"mtx_species_sample_distance_in_PCA.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
  #save the corresponding eigenvalue ratios
write.table(x=df_sample_pc_eigenvalues,file = paste0(output_workspace,"Table_df_sample_pc_eigenvalues.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
  #save the corresponding PCA projection
write.table(x=df_sample_pca,file = paste0(output_workspace,"Table_df_sample_pca.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
  #save the corresponding PC loadings 
write.table(x=sample_pca$rotation,file = paste0(output_workspace,"mtx_species_sample_pca_loadings.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_sorted_pc1_and_pc2_loadings,file = paste0(output_workspace,"Table_df_sorted_pc1_and_pc2_loadings.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

#Multiple linear regression sp_samples_pca_distance ~ species_phylo_distance + is_from_same_condition + is_from_same_sample_replicate + is_from_same_mixture
df_species_sample_pca_dist_and_features <- as.data.frame(as_tibble(mtx_species_sample_distance_in_PCA, rownames = "species_sample_1")) %>% 
  #convert to a "long" format
  pivot_longer(cols=-1, names_to = "species_sample_2", values_to = "sp_samples_pca_distance")
  #Add features
    #add corresponding species pair phylogenetic distance
v_sp1_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = v_species_sample_name_to_species_smpl_lbl[df_species_sample_pca_dist_and_features$species_sample_1]))
v_sp2_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = v_species_sample_name_to_species_smpl_lbl[df_species_sample_pca_dist_and_features$species_sample_2]))
df_species_sample_pca_dist_and_features$species_phylo_distance <- vapply(X = 1:length(v_sp1_in_pair),FUN = function(the_i) mtx_phylo_dist_the_21_sps[v_sp1_in_pair[the_i],v_sp2_in_pair[the_i]],FUN.VALUE = 0.0)
    #add corresponding is_from_same_condition
df_species_sample_pca_dist_and_features$sp_sample1_condition <- substr(x = df_species_sample_pca_dist_and_features$species_sample_1,start = 4,stop = 6)
df_species_sample_pca_dist_and_features$sp_sample2_condition <- substr(x = df_species_sample_pca_dist_and_features$species_sample_2,start = 4,stop = 6)
df_species_sample_pca_dist_and_features$is_from_same_condition <- df_species_sample_pca_dist_and_features$sp_sample1_condition == df_species_sample_pca_dist_and_features$sp_sample2_condition
    #add corresponding is_from_same_sample_replicate
df_species_sample_pca_dist_and_features$sp_sample1_replicate <- gsub(pattern = "RNA_1_",replacement = "",x = gsub(pattern = "RNA_2_",replacement = "",x = df_species_sample_pca_dist_and_features$species_sample_1))
df_species_sample_pca_dist_and_features$sp_sample2_replicate <- gsub(pattern = "RNA_1_",replacement = "",x = gsub(pattern = "RNA_2_",replacement = "",x = df_species_sample_pca_dist_and_features$species_sample_2))
df_species_sample_pca_dist_and_features$sp_sample1_replicate <- stringr::str_remove_all(df_species_sample_pca_dist_and_features$sp_sample1_replicate, "_S\\d+")
df_species_sample_pca_dist_and_features$sp_sample2_replicate <- stringr::str_remove_all(df_species_sample_pca_dist_and_features$sp_sample2_replicate, "_S\\d+")
df_species_sample_pca_dist_and_features$is_from_same_sample_replicate <- df_species_sample_pca_dist_and_features$sp_sample1_replicate == df_species_sample_pca_dist_and_features$sp_sample2_replicate
    #add corresponding is_from_same_mixture
df_species_sample_pca_dist_and_features$sp_sample1_mixture <- substr(x = df_species_sample_pca_dist_and_features$species_sample_1,start = 1,stop = 2)
df_species_sample_pca_dist_and_features$sp_sample2_mixture <- substr(x = df_species_sample_pca_dist_and_features$species_sample_2,start = 1,stop = 2)
df_species_sample_pca_dist_and_features$is_from_same_mixture <- df_species_sample_pca_dist_and_features$sp_sample1_mixture == df_species_sample_pca_dist_and_features$sp_sample2_mixture
    #remove sample pairs with a sample compare to itself
df_species_sample_pca_dist_and_features <- subset(df_species_sample_pca_dist_and_features,(species_sample_1!=species_sample_2))

  #run multiple linear regression PCA distance ~ species_phylo_distance + is_from_same_condition + is_from_same_mixture 
lm_sp_sample_pca_dist <- lm(data = df_species_sample_pca_dist_and_features, formula = sp_samples_pca_distance ~ species_phylo_distance + is_from_same_condition + is_from_same_mixture)
summary(lm_sp_sample_pca_dist) #There is no batch effect as the similarity of mixtures does not signficantly affect the expression PCA distance between sample pairs
coefficients()
  #Add the sample pairs expression data R^2 to the dataframe
mtx_all_66samples_expr_corr_R2 <- as.matrix(cor(t(mtx_log2_1p_sp_sample_FPKM))^2)
df_species_sample_pca_dist_and_features$sample_pair_expression_R2 <- vapply(X = 1:nrow(df_species_sample_pca_dist_and_features),FUN = function(the_i) mtx_all_66samples_expr_corr_R2[df_species_sample_pca_dist_and_features$species_sample_1[the_i],df_species_sample_pca_dist_and_features$species_sample_2[the_i]],FUN.VALUE = c(0.0))
    #Visualize R2 matrix
corrplot::corrplot(mtx_all_66samples_expr_corr_R2, type = "upper", method = "color", tl.pos = "n")
row_dist <- vegdist(mtx_all_66samples_expr_corr_R2, method = "euclidean")
col_dist <- vegdist(t(mtx_all_66samples_expr_corr_R2), method = "euclidean")
row_clust <- hclust(row_dist, method = "ward.D2") #Use ward.D2 for Ward's method
col_clust <- hclust(col_dist, method = "ward.D2")
png(filename = paste0(output_workspace,"Heatmap_all_66samples_expr_corr_R2.png"),width = 25,height = 25,units = "cm",res = 300)
heatmap.2(mtx_all_66samples_expr_corr_R2,
                                             distfun = function(x) vegdist(x, method = "euclidean"), #Specify Bray-Curtis for heatmap.2
                                             hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
                                             Rowv = as.dendrogram(row_clust),
                                             Colv = as.dendrogram(col_clust),
                                             dendrogram = "none",
                                             trace = "none",
                                             main = NULL,
                                             xlab = "Sample 1",
                                             ylab="Sample 2",cexRow=0.4,cexCol=0.4,
                                             key.par = list(cex = 0.45),
                                             key.title = NA,
                                             key.xlab = "log10(pathway size+1)",
                                             key.ylab = "Count",
                                             col= brewer.pal(n = 9, name = "YlOrRd"),
                                             margins = c(5.5, 15.5),
                                             lwid = c(1.2, 8),
                                             lhei = c(0.5, 2))
dev.off()
    #run multiple linear regression sample_pair_expression_R2 ~ species_phylo_distance + is_from_same_condition + is_from_same_mixture
lm_sample_pairs_expr_R2_vs_sample_features <- lm(data = df_species_sample_pca_dist_and_features, formula = sample_pair_expression_R2 ~ species_phylo_distance + is_from_same_condition  + is_from_same_mixture)
summary(lm_sample_pairs_expr_R2_vs_sample_features) #There is no batch effect as the similarity of mixtures does not signficantly affect the expression correlation (R2) between sample pairs

  #Compare the expression R2 of the YPD sample pairs from the same mixture for pairs from the same species vs pairs from different species
    #remove samples from other conditions than YPD and sample pairs with different conditions 
df_species_sample_pca_dist_and_features_YPD <- subset(df_species_sample_pca_dist_and_features,(sp_sample1_condition=="YPD")&(sp_sample2_condition=="YPD"))
    #now that we filtered for the right samples, define whether or not samples are from the same species
df_species_sample_pca_dist_and_features_YPD$is_from_same_species <- df_species_sample_pca_dist_and_features_YPD$species_phylo_distance == 0
    #compare R2 of sample pairs from YPD, the same species and the same mixtures (replicates) TO R2 of sample pairs from the same
ggboxplot(df_species_sample_pca_dist_and_features_YPD, x = "is_from_same_species", y = "sample_pair_expression_R2") +
  stat_compare_means(
    method = "wilcox.test", 
    method.args = list(alternative = "greater"),label = "p.format",
    label.y = 1  #Adjust position of the p-value label
  ) + xlab("Are the samples technical replicates?") + ylab(expression("Sample pairs expression " * R^2)) + scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none")
ggsave(filename = "Sample_pairs_expr_Rsq_for_replicates_vs_nonreplicates_in_YPD.png", path=output_workspace, width = 15/2.54, height = 10/2.54)
ggsave(filename = "Sample_pairs_expr_Rsq_for_replicates_vs_nonreplicates_in_YPD.svg", path=output_workspace, width = 15/2.54, height = 10/2.54,device=svg)

#Create a PCA of all species aggregated FPKM
  #get species aggregate matrices (best aggregate function is mean)
mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YPD"), species_smpl_lbl ~ Orthogroup, value.var = "log2_1p_FPKM",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD,file = paste0(output_workspace,"mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YNB"), species_smpl_lbl ~ Orthogroup, value.var = "log2_1p_FPKM",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB,file = paste0(output_workspace,"mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

  #Perform the PCA for both conditions
species_aggregate_gene_expr_pca_YPD <- prcomp(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD)
species_aggregate_gene_expr_pca_YNB <- prcomp(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB)
  #Visualize PCA
species_aggregate_gene_expr_pc_eigenvalues_YPD <- species_aggregate_gene_expr_pca_YPD$sdev^2
species_aggregate_gene_expr_pc_eigenvalues_YNB <- species_aggregate_gene_expr_pca_YNB$sdev^2

    #YPD
df_species_aggregate_gene_expr_pca_YPD <- as.data.frame(as_tibble(species_aggregate_gene_expr_pca_YPD$x, rownames = "species"))
df_species_aggregate_gene_expr_pca_YPD$species_initials <- vapply(X = df_species_aggregate_gene_expr_pca_YPD$species,FUN = get_initials,FUN.VALUE = "")
df_species_aggregate_gene_expr_pca_YPD$mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_YPD$species]
      #Eigenvalues
species_aggregate_gene_expr_pc_eigenvalues_YPD <- species_aggregate_gene_expr_pca_YPD$sdev^2
df_species_aggregate_gene_expr_pc_eigenvalues_YPD <- tibble(PC = factor(1:length(species_aggregate_gene_expr_pc_eigenvalues_YPD)), 
                                   variance = species_aggregate_gene_expr_pc_eigenvalues_YPD) %>% 
        #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
        #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YPD <- as.integer(subset(df_species_aggregate_gene_expr_pc_eigenvalues_YPD,explained_variance_cum>=0.99)$PC)[1]
        #Visualize species aggregate gene expression PCA eigen values
df_species_aggregate_gene_expr_pc_eigenvalues_YPD %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = explained_variance)) +
  geom_line(aes(y = explained_variance_cum, group = 1)) + 
  geom_point(aes(y = explained_variance_cum)) +
  geom_hline(yintercept = 0.99,col="red",lty=2)+
  geom_vline(xintercept = nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YPD,col="red",lty=2)+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  labs(x = "Principal component", y = "Fraction variance explained")
ggsave(filename = "PCA_species_gene_expr_aggregate_Cumul_Explained_variance_YPD.png", path=output_workspace, width = 25/2.54, height = 15/2.54)

        #Visualize species aggregate gene expression PCA
df_species_aggregate_gene_expr_pca_YPD %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(size = 2, stroke = 2) + 
  ggrepel::geom_text_repel(aes(label = mixture)) #geom_text
ggsave(filename = "PCA_species_gene_expr_aggregate_YPD.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

        #Visualize top 10 ORTHOGROUPS (or "gene") with highest loadings in species aggregate gene expression PCA's pc1 and pc2
species_aggregate_gene_expr_pc_loadings_YPD <- species_aggregate_gene_expr_pca_YPD$rotation %>% 
  as_tibble(rownames = "gene") %>%
  as.data.frame()

v_species_aggregate_gene_expr_pc_eigenvalues_YPD <- species_aggregate_gene_expr_pc_eigenvalues_YPD
names(v_species_aggregate_gene_expr_pc_eigenvalues_YPD) <- paste0("PC",1:length(species_aggregate_gene_expr_pc_eigenvalues_YPD))

#as pc1/pc2 loading dataframe
df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD <- species_aggregate_gene_expr_pc_loadings_YPD %>% 
  #select only the PCs we are interested in
  dplyr::select(gene, PC1, PC2) %>%
  #convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  #for each PC
  dplyr::group_by(PC) %>%
  as.data.frame()

df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$the_princ_comp_explained_var <- v_species_aggregate_gene_expr_pc_eigenvalues_YPD[df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$PC]
df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$normalized_loading <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$loading*df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$the_princ_comp_explained_var

df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD %>% 
  #arrange by descending order of loading
  dplyr::arrange(desc(abs(normalized_loading))) 

species_aggregate_gene_expr_top10_gene_orthogroup_YPD <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD$gene[1:10]

species_aggregate_gene_expr_top10_gene_orthogroup_YPD
#Functions species_aggregate_gene_expr_top10_gene_orthogroup_YPD:
species_aggregate_gene_expr_top10_gene_orthogroup_functions_YPD <- c("Cell wall oganization protein Zeo1", "NCW1 protein of unknown function (the orthogroup is Haniespora-specific although another orthogroup shares the same functional annotation in other species too and is the top10th orthgroup)", "unknown", "Cell wall stabilizing protein CWP2 involved in the resistance to low pH", "unknown but specific to Hanseniaspora", "STF2 ATPase involved in stress resistance (oxidative and dessication)", "posttranslational protein targeting EGD2", "unknown but specific to Hanseniaspora", "Mitochondrial ATPase inhibitor, IATP", "NCW1 protein of unknown function")
species_aggregate_gene_expr_top10_gene_orthogroup_functions_YPD

      #visualize sample PC1/PC2 top 10 orthogroup pc loading
species_aggregate_gene_expr_top_loadings_YPD <- species_aggregate_gene_expr_pc_loadings_YPD %>% 
  dplyr::filter(gene %in% species_aggregate_gene_expr_top10_gene_orthogroup_YPD)

ggplot(data = species_aggregate_gene_expr_top_loadings_YPD) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               col = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
                           nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))+
  labs(x = "PC1", y = "PC2")
ggsave(filename = "Top10_pc1_pc2_normalized_loadings_PCA_species_gene_expr_aggregate_YPD.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

  #YNB
df_species_aggregate_gene_expr_pca_YNB <- as.data.frame(as_tibble(species_aggregate_gene_expr_pca_YNB$x, rownames = "species"))
df_species_aggregate_gene_expr_pca_YNB$species_initials <- vapply(X = df_species_aggregate_gene_expr_pca_YNB$species,FUN = get_initials,FUN.VALUE = "")
df_species_aggregate_gene_expr_pca_YNB$mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_YNB$species]
    #Eigenvalues
species_aggregate_gene_expr_pc_eigenvalues_YNB <- species_aggregate_gene_expr_pca_YNB$sdev^2
df_species_aggregate_gene_expr_pc_eigenvalues_YNB <- tibble(PC = factor(1:length(species_aggregate_gene_expr_pc_eigenvalues_YNB)), 
  variance = species_aggregate_gene_expr_pc_eigenvalues_YNB) %>% 
    #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
    #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YNB <- as.integer(subset(df_species_aggregate_gene_expr_pc_eigenvalues_YNB,explained_variance_cum>=0.99)$PC)[1]
    #Visualize species aggregate gene expression PCA eigen values
df_species_aggregate_gene_expr_pc_eigenvalues_YNB %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = explained_variance)) +
  geom_line(aes(y = explained_variance_cum, group = 1)) + 
  geom_point(aes(y = explained_variance_cum)) +
  geom_hline(yintercept = 0.99,col="red",lty=2)+
  geom_vline(xintercept = nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YNB,col="red",lty=2)+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  labs(x = "Principal component", y = "Fraction variance explained")
ggsave(filename = "PCA_species_gene_expr_aggregate_Cumul_Explained_variance_YNB.png", path=output_workspace, width = 25/2.54, height = 15/2.54)

    #Visualize species aggregate gene expression PCA
df_species_aggregate_gene_expr_pca_YNB %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(size = 2, stroke = 2) + 
  ggrepel::geom_text_repel(aes(label = mixture)) #geom_text
ggsave(filename = "PCA_species_gene_expr_aggregate_YNB.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

    #Visualize top 10 ORTHOGROUPS (or "gene") with highest loadings in species aggregate gene expression PCA's pc1 and pc2
species_aggregate_gene_expr_pc_loadings_YNB <- species_aggregate_gene_expr_pca_YNB$rotation %>% 
  as_tibble(rownames = "gene") %>%
  as.data.frame()

v_species_aggregate_gene_expr_pc_eigenvalues_YNB <- species_aggregate_gene_expr_pc_eigenvalues_YNB
names(v_species_aggregate_gene_expr_pc_eigenvalues_YNB) <- paste0("PC",1:length(species_aggregate_gene_expr_pc_eigenvalues_YNB))

    #as pc1/pc2 loading dataframe
df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB <- species_aggregate_gene_expr_pc_loadings_YNB %>% 
      #select only the PCs we are interested in
  dplyr::select(gene, PC1, PC2) %>%
      #convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
      #for each PC
  dplyr::group_by(PC) %>%
  as.data.frame()

df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$the_princ_comp_explained_var <- v_species_aggregate_gene_expr_pc_eigenvalues_YNB[df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$PC]
df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$normalized_loading <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$loading*df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$the_princ_comp_explained_var

df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB %>% 
      #arrange by descending order of loading
  dplyr::arrange(desc(abs(normalized_loading))) 

species_aggregate_gene_expr_top10_gene_orthogroup_YNB <- df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB$gene[1:10]

species_aggregate_gene_expr_top10_gene_orthogroup_YNB
#Functions species_aggregate_gene_expr_top10_gene_orthogroup_YNB:
species_aggregate_gene_expr_top10_gene_orthogroup_functions_YNB <- c("Cell wall organization mannoprotein CWP2",
                                                                     "Cell membrane channel/ Gap junction alpha-8 protein GJA8",
                                                                     "unknown",
                                                                     "NCW1 protein of unknown function",
                                                                     "unknown",
                                                                     "putative phosphatase",
                                                                     "unknown",
                                                                     "unknown",
                                                                     "unknown",
                                                                     "DPI8 Delta-Psi dependent mitochondrial Import protein of 8 kDa")
species_aggregate_gene_expr_top10_gene_orthogroup_functions_YNB

    #visualize sample PC1/PC2 top 10 orthogroup pc loading
species_aggregate_gene_expr_top_loadings_YNB <- species_aggregate_gene_expr_pc_loadings_YNB %>% 
  dplyr::filter(gene %in% species_aggregate_gene_expr_top10_gene_orthogroup_YNB)

ggplot(data = species_aggregate_gene_expr_top_loadings_YNB) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               col = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
                           nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))+
  labs(x = "PC1", y = "PC2")
ggsave(filename = "Top10_pc1_pc2_normalized_loadings_PCA_species_gene_expr_aggregate_YNB.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

#Get species PCAs pairwise distances
  #max
mtx_species_aggregate_gene_expr_pca_distance_YPD <- as.matrix(vegan::vegdist(species_aggregate_gene_expr_pca_YPD$x, method = "euclidean"))
mtx_species_aggregate_gene_expr_pca_distance_YNB <- as.matrix(vegan::vegdist(species_aggregate_gene_expr_pca_YNB$x, method = "euclidean"))

#create useful conversion vectors
v_species_smpl_lbl_to_mixture <- v_species_to_mixture
names(v_species_smpl_lbl_to_mixture) <- v_species_to_simplified_species_lbl[names(v_species_smpl_lbl_to_mixture)]

#Simple linear regression species_aggregate_gene_expr_pca_distance ~ species_phylo_distance
  #FOR YPD
df_species_aggregate_gene_expr_pca_dist_and_features_YPD <- as.data.frame(as_tibble(mtx_species_aggregate_gene_expr_pca_distance_YPD, rownames = "species1")) %>% 
    #convert to a "long" format
  pivot_longer(cols=-1, names_to = "species2", values_to = "species_aggregate_gene_expr_pca_distance")
    #Add features
    #add corresponding species pair phylogenetic distance
v_sp1_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species1))
v_sp2_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species2))
df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species_phylo_distance <- vapply(X = 1:length(v_sp1_in_pair),FUN = function(the_i) mtx_phylo_dist_the_21_sps[v_sp1_in_pair[the_i],v_sp2_in_pair[the_i]],FUN.VALUE = 0.0)
    #add corresponding is_from_same_mixture
df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species1_mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species1]
df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species2_mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species2]
df_species_aggregate_gene_expr_pca_dist_and_features_YPD$is_from_same_mixture <- df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species1_mixture == df_species_aggregate_gene_expr_pca_dist_and_features_YPD$species2_mixture
    #run multiple linear regression
lm_species_aggregate_gene_expr_pca_dist_YPD <- lm(data = df_species_aggregate_gene_expr_pca_dist_and_features_YPD, formula = species_aggregate_gene_expr_pca_distance ~ species_phylo_distance + is_from_same_mixture)
summary(lm_species_aggregate_gene_expr_pca_dist_YPD)
    #run simple linear regression
simple_lm_species_aggregate_gene_expr_pca_dist_YPD <- lm(data = df_species_aggregate_gene_expr_pca_dist_and_features_YPD, formula = species_aggregate_gene_expr_pca_distance ~ species_phylo_distance)
summary(simple_lm_species_aggregate_gene_expr_pca_dist_YPD)

  #FOR YNB
df_species_aggregate_gene_expr_pca_dist_and_features_YNB <- as.data.frame(as_tibble(mtx_species_aggregate_gene_expr_pca_distance_YNB, rownames = "species1")) %>% 
    #convert to a "long" format
  pivot_longer(cols=-1, names_to = "species2", values_to = "species_aggregate_gene_expr_pca_distance")
    #Add features
    #add corresponding species pair phylogenetic distance
v_sp1_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species1))
v_sp2_in_pair <- unname(gsub(pattern = " ",replacement = "_",x = df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species2))
df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species_phylo_distance <- vapply(X = 1:length(v_sp1_in_pair),FUN = function(the_i) mtx_phylo_dist_the_21_sps[v_sp1_in_pair[the_i],v_sp2_in_pair[the_i]],FUN.VALUE = 0.0)
    #add corresponding is_from_same_mixture
df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species1_mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species1]
df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species2_mixture <- v_species_smpl_lbl_to_mixture[df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species2]
df_species_aggregate_gene_expr_pca_dist_and_features_YNB$is_from_same_mixture <- df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species1_mixture == df_species_aggregate_gene_expr_pca_dist_and_features_YNB$species2_mixture
    #run multiple linear regression
lm_species_aggregate_gene_expr_pca_dist_YNB <- lm(data = df_species_aggregate_gene_expr_pca_dist_and_features_YNB, formula = species_aggregate_gene_expr_pca_distance ~ species_phylo_distance + is_from_same_mixture)
summary(lm_species_aggregate_gene_expr_pca_dist_YNB)
    #run simple linear regression
simple_lm_species_aggregate_gene_expr_pca_dist_YNB <- lm(data = df_species_aggregate_gene_expr_pca_dist_and_features_YNB, formula = species_aggregate_gene_expr_pca_distance ~ species_phylo_distance)
summary(simple_lm_species_aggregate_gene_expr_pca_dist_YNB)

    #save the PCA stats as tsv files
      #save pca matrix
write.table(x=species_aggregate_gene_expr_pca_YPD$x,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_pca_YNB$x,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
      #save pca disance matrix
write.table(x=mtx_species_aggregate_gene_expr_pca_distance_YPD,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_distance_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=mtx_species_aggregate_gene_expr_pca_distance_YNB,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_distance_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
      #save the corresponding eigenvalue ratios
write.table(x=df_species_aggregate_gene_expr_pc_eigenvalues_YPD,file = paste0(output_workspace,"Table_df_species_aggregate_gene_expr_pc_eigenvalues_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_species_aggregate_gene_expr_pc_eigenvalues_YNB,file = paste0(output_workspace,"Table_df_species_aggregate_gene_expr_pc_eigenvalues_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
      #save the corresponding PCA projection
write.table(x=df_species_aggregate_gene_expr_pca_YPD,file = paste0(output_workspace,"Table_df_species_aggregate_gene_expr_pca_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_species_aggregate_gene_expr_pca_YNB,file = paste0(output_workspace,"Table_df_species_aggregate_gene_expr_pca_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
      #save the corresponding PC loadings 
write.table(x=species_aggregate_gene_expr_pca_YPD$rotation,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_loadings_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YPD,file = paste0(output_workspace,"Table_df_sorted_pc1_and_pc2_loadings_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_pca_YNB$rotation,file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_loadings_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_sorted_species_aggregate_gene_expr_pc1_and_pc2_loadings_YNB,file = paste0(output_workspace,"Table_df_sorted_pc1_and_pc2_loadings_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

#Pairwise expression distance matrices (Euclidean)
mtx_pairwise_expression_dist_YPD <- as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD, method = "euclidean"))
mtx_pairwise_expression_dist_YNB <- as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB, method = "euclidean"))

#Pairwise expression similarity/correlation matrices (R^2)
mtx_pairwise_expression_similarity_YPD <- as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD))^2)
mtx_gene_content_with_NA_for_absence_in_species_genomes_YPD <- mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD
mtx_gene_content_with_NA_for_absence_in_species_genomes_YPD[mtx_gene_content_with_NA_for_absence_in_species_genomes_YPD==0] <- NA
mtx_pairwise_expression_similarity_based_on_present_genes_YPD <- as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD*mtx_gene_content_with_NA_for_absence_in_species_genomes_YPD),use = "pairwise.complete.obs")^2)
mtx_pairwise_expression_similarity_YNB <- as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB))^2)
mtx_gene_content_with_NA_for_absence_in_species_genomes_YNB <- mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB
mtx_gene_content_with_NA_for_absence_in_species_genomes_YNB[mtx_gene_content_with_NA_for_absence_in_species_genomes_YNB==0] <- NA
mtx_pairwise_expression_similarity_based_on_present_genes_YNB <- as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB*mtx_gene_content_with_NA_for_absence_in_species_genomes_YNB),use = "pairwise.complete.obs")^2)

  #Cast df_yeasts_pantranscriptome to genomic presence/absence matrix
mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YPD"), species_smpl_lbl ~ Orthogroup, value.var = "int_is_detected_in_genome",fun.aggregate = function(x) max(x,na.rm=T),fill = 0)
write.table(x=mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD,file = paste0(output_workspace,"mtx_gene_pres_abs_in_species_genomes_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YNB"), species_smpl_lbl ~ Orthogroup, value.var = "int_is_detected_in_genome",fun.aggregate = function(x) max(x,na.rm=T),fill = 0)
write.table(x=mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB,file = paste0(output_workspace,"mtx_gene_pres_abs_in_species_genomes_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
mtx_gene_INTEGER_pres_abs_in_species_genomes <- reshape2::acast(df_yeasts_pantranscriptome, species_smpl_lbl ~ Orthogroup, value.var = "int_is_detected_in_genome",fun.aggregate = function(x) max(x,na.rm=T),fill = 0)
write.table(x=mtx_gene_INTEGER_pres_abs_in_species_genomes,file = paste0(output_workspace,"mtx_gene_pres_abs_in_species_genomes.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
proportion_core_genes_in_the_21_species_pangenome <- sum(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))==1)/nrow(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))
proportion_core_genes_in_the_21_species_pangenome
proportion_accessory_genes_in_the_21_species_pangenome <- 1-sum(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))==1)/nrow(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))
proportion_accessory_genes_in_the_21_species_pangenome
#sort((rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes)))[(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes)))<1],dec=T)
mtx_gene_BOOL_pres_abs_in_species_genomes_YPD <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YPD"), species_smpl_lbl ~ Orthogroup, value.var = "int_is_detected_in_genome",fun.aggregate = function(x) max(x,na.rm=T),fill = 0)
mtx_gene_BOOL_pres_abs_in_species_genomes_YNB <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YNB"), species_smpl_lbl ~ Orthogroup, value.var = "int_is_detected_in_genome",fun.aggregate = function(x) max(x,na.rm=T),fill = 0)
mtx_BOOL_is_gene_expressed_in_species_in_YPD <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YPD"), species_smpl_lbl ~ Orthogroup, value.var = "bool_is_expressed",fun.aggregate = function(x) as.logical(max(x,na.rm=T)),fill = FALSE)
mtx_BOOL_is_gene_expressed_in_species_in_YNB <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YNB"), species_smpl_lbl ~ Orthogroup, value.var = "bool_is_expressed",fun.aggregate = function(x) as.logical(max(x,na.rm=T)),fill = FALSE)

#save gene/orthogroup presence/absence matrix transpose across species
write.table(x=t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD),file = paste0(output_workspace,"mtx_TRANSPOSE_gene_pres_abs_in_species_genomes_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB),file = paste0(output_workspace,"mtx_TRANSPOSE_gene_pres_abs_in_species_genomes_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

#PCA on species gene expression aggregate matrix transpose (so that genes/orthogroups can be used as the rows/objects of the varpart)
  #Perform the PCA for both conditions
species_aggregate_gene_expr_transpose_pca_YPD <- prcomp(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD))
species_aggregate_gene_expr_transpose_pca_YNB <- prcomp(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB))
  #Visualize PCA
species_aggregate_gene_expr_transpose_pca_pc_eigenvalues_YPD <- species_aggregate_gene_expr_transpose_pca_YPD$sdev^2
species_aggregate_gene_expr_transpose_pca_pc_eigenvalues_YNB <- species_aggregate_gene_expr_transpose_pca_YNB$sdev^2

    #YPD
df_species_aggregate_gene_expr_transpose_pca_YPD <- as.data.frame(as_tibble(species_aggregate_gene_expr_transpose_pca_YPD$x, rownames = "Genes_orthogroup"))
      #Eigenvalues
species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD <- species_aggregate_gene_expr_transpose_pca_YPD$sdev^2
df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD <- tibble(PC = factor(1:length(species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD)), 
                                                            variance = species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD) %>% 
        #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
        #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

#View(species_aggregate_gene_expr_transpose_pca_YPD$rotation)

nb_pcs_recapitulating_species_aggregate_gene_expr_transpose_pca_YPD <- as.integer(subset(df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD,explained_variance_cum>=0.99)$PC)[1]

    #YNB
df_species_aggregate_gene_expr_transpose_pca_YNB <- as.data.frame(as_tibble(species_aggregate_gene_expr_transpose_pca_YNB$x, rownames = "Genes_orthogroup"))
      #Eigenvalues
species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB <- species_aggregate_gene_expr_transpose_pca_YNB$sdev^2
df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB <- tibble(PC = factor(1:length(species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB)), 
                                                                      variance = species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB) %>% 
        #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
        #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

#View(species_aggregate_gene_expr_transpose_pca_YNB$rotation)

nb_pcs_recapitulating_species_aggregate_gene_expr_transpose_pca_YNB <- as.integer(subset(df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB,explained_variance_cum>=0.99)$PC)[1]

#save species gene expression aggregate matrix transpose PCA
write.table(x=df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD,file = paste0(output_workspace,"Table_df_TRANSPOSE_species_aggregate_gene_expr_pc_eigenvalues_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB,file = paste0(output_workspace,"Table_df_TRANSPOSE_species_aggregate_gene_expr_pc_eigenvalues_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_transpose_pca_YPD$x,file = paste0(output_workspace,"mtx_TRANSPOSE_species_aggregate_gene_expr_pca_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_transpose_pca_YNB$x,file = paste0(output_workspace,"mtx_TRANSPOSE_species_aggregate_gene_expr_pca_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_transpose_pca_YPD$x[,1:nb_pcs_recapitulating_species_aggregate_gene_expr_transpose_pca_YPD],file = paste0(output_workspace,"mtx_TRANSPOSE_species_aggregate_gene_expr_pca_RETAINED_PCs_YPD.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_transpose_pca_YNB$x[,1:nb_pcs_recapitulating_species_aggregate_gene_expr_transpose_pca_YNB],file = paste0(output_workspace,"mtx_TRANSPOSE_species_aggregate_gene_expr_pca_RETAINED_PCs_YNB.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

#Get the %variance of expression explained by gene presence/absence (pangenomic component)
#  #Based on pairwise distances
#mtx_dist_species_aggregate_expr_YPD <- as.matrix(1/mtx_pairwise_expression_similarity_YPD)#as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD,method = "euclidean"))
#mtx_dist_species_pres_abs_YPD <- as.matrix(vegan::vegdist(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),],method = "gower"))
#v_dist_species_aggregate_expr_YPD <- NULL
#v_dist_species_pres_abs_YPD <- NULL
#for (the_i in 1:nrow(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD)){
#  for (the_j in 1:the_i){
#    v_dist_species_aggregate_expr_YPD <- c(v_dist_species_aggregate_expr_YPD, mtx_dist_species_aggregate_expr_YPD[the_i,the_j])
#    v_dist_species_pres_abs_YPD <- c(v_dist_species_pres_abs_YPD,mtx_dist_species_pres_abs_YPD[the_i,the_j])
#  }
#}
#
#mtx_dist_species_aggregate_expr_YNB <- as.matrix(1/mtx_pairwise_expression_similarity_YNB)#as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB,method = "euclidean"))
#mtx_dist_species_pres_abs_YNB <- as.matrix(vegan::vegdist(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB[rownames(mtx_pairwise_expression_similarity_YNB),],method = "gower"))
#v_dist_species_aggregate_expr_YNB <- NULL
#v_dist_species_pres_abs_YNB <- NULL
#for (the_i in 1:nrow(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB)){
#  for (the_j in 1:the_i){
#    v_dist_species_aggregate_expr_YNB <- c(v_dist_species_aggregate_expr_YNB, mtx_dist_species_aggregate_expr_YNB[the_i,the_j])
#    v_dist_species_pres_abs_YNB <- c(v_dist_species_pres_abs_YNB,mtx_dist_species_pres_abs_YNB[the_i,the_j])
#  }
#}
#    #YPD
#summary(lm(v_dist_species_aggregate_expr_YPD~(v_dist_species_pres_abs_YPD) ))$r.squared
#    #YNB
#summary(lm(v_dist_species_aggregate_expr_YNB~(v_dist_species_pres_abs_YNB) ))$r.squared

  #Based on multiple linear regression (MLR)
    #YPD
v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD <- df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YPD$explained_variance
v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD <- rep(NA,length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD))
for (the_pc_id in 1:length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD)){
  v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD[the_pc_id] <- summary(lm(species_aggregate_gene_expr_transpose_pca_YPD$x[,the_pc_id]~t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD)))$r.squared
}
var_ratio_of_expression_explained_by_presence_absence_YPD_MLR <- sum(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD*v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD)
var_ratio_of_expression_explained_by_presence_absence_YPD_MLR

      #bootstrap var ratio of expression_explained_by_presence_absence_YPD (999 permutations of 100 genes)
nb_subsamples <- 999
v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YPD_MLR <- NULL
for (the_current_permutation_id in 1:nb_subsamples){
  v_current_permutation_gene_indices <- sample(x = 1:nrow(species_aggregate_gene_expr_transpose_pca_YPD$x),size = 100,replace = F)
  v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD <- rep(NA,length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD))
  for (the_pc_id in 1:length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD)){
    v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD[the_pc_id] <- summary(lm(species_aggregate_gene_expr_transpose_pca_YPD$x[v_current_permutation_gene_indices,the_pc_id]~t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD)[v_current_permutation_gene_indices,]))$r.squared
  }
  v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YPD_MLR <- c(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YPD_MLR,sum(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YPD*v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YPD))
  if (the_current_permutation_id%%100==0){
    print(the_current_permutation_id)
  }
}
summary(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YPD_MLR);sd(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YPD_MLR)

    #YNB
v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB <- df_species_aggregate_gene_expr_transpose_pc_eigenvalues_YNB$explained_variance
v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB <- rep(NA,length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB))
for (the_pc_id in 1:length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB)){
  v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB[the_pc_id] <- summary(lm(species_aggregate_gene_expr_transpose_pca_YNB$x[,the_pc_id]~t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB)))$r.squared
}
var_ratio_of_expression_explained_by_presence_absence_YNB_MLR <- sum(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB*v_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB)
var_ratio_of_expression_explained_by_presence_absence_YNB_MLR
      #bootstrap var ratio of expression_explained_by_presence_absence_YNB
nb_subsamples <- 999
v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YNB_MLR <- NULL
for (the_current_permutation_id in 1:nb_subsamples){
  v_current_permutation_gene_indices <- sample(x = 1:nrow(species_aggregate_gene_expr_transpose_pca_YNB$x),size = 100,replace = F)
  v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB <- rep(NA,length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB))
  for (the_pc_id in 1:length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB)){
    v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB[the_pc_id] <- summary(lm(species_aggregate_gene_expr_transpose_pca_YNB$x[v_current_permutation_gene_indices,the_pc_id]~t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB)[v_current_permutation_gene_indices,]))$r.squared
  }
  v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YNB_MLR <- c(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YNB_MLR,sum(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB*v_BOOTSTRAP_species_aggregate_gene_expr_Transpose_pc_vs_Pres_abs_R2_YNB))
  if (the_current_permutation_id%%100==0){
    print(the_current_permutation_id)
  }
}
summary(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YNB_MLR);sd(v_BOOTSTRAP_var_ratio_of_expression_explained_by_presence_absence_YNB_MLR)

  #Based on linear mixed model
#    #YPD
#var_ratio_of_expression_explained_by_presence_absence_YPD_LMM <- sum(as.numeric(read.csv2(file = paste0(output_workspace,"Table_weighted_Expression_PCs_R2_corr_with_The_Pres_Abs_component_YPD.tsv"),sep = "\t",header = F,stringsAsFactors = FALSE)[,3]))/0.9
#    #YNB
#var_ratio_of_expression_explained_by_presence_absence_YNB_LMM <- sum(as.numeric(read.csv2(file = paste0(output_workspace,"Table_weighted_Expression_PCs_R2_corr_with_The_Pres_Abs_component_YNB.tsv"),sep = "\t",header = F,stringsAsFactors = FALSE)[,3]))/0.9

#S. cerevisiae within-species %variance of expression explained by gene presence/absence
  #Import Scer strains gene presence/absence data
org_mtx_Scer_pres_abs <- read.csv2(file = paste0(output_workspace,"Sc_genesMatrix_PresenceAbsence.tab"),sep = "\t",row.names = 1,header=T,stringsAsFactors = FALSE)

#import Sc expression data and define matchable gene names
df_Scer_expression <- read.csv2(file = paste0(output_workspace,"Sc_strains_Expression_data.tab"),sep = ",",header = T,stringsAsFactors = FALSE)
df_Scer_expression$tpm <- as.numeric(df_Scer_expression$tpm)
df_Scer_expression$matchable_gene_name <- vapply(X = df_Scer_expression$ORF,FUN = function(the_nm) unlist( strsplit(x = the_nm,split="[.]"))[1],FUN.VALUE = "")

#create the gene content matrix and fix the gene names 
mtx_Scer_pres_abs <- org_mtx_Scer_pres_abs
v_gene_nm_start_Y <- gsub(pattern = ".",replacement = "",x = str_extract(string = colnames(mtx_Scer_pres_abs),pattern = "\\.Y.{6}\\.[A-Z]|\\.Y.{6}"),fixed = T)
v_gene_nm_start_X <- gsub(pattern = ".",replacement = "",x = str_extract(string = colnames(mtx_Scer_pres_abs),pattern = "X.*?\\."),fixed = T)

for (j in 1:ncol(mtx_Scer_pres_abs)){
  current_gene_nm <- colnames(mtx_Scer_pres_abs)[j]
  current_gene_nm_start_Y <- v_gene_nm_start_Y[j]#gsub(pattern = ".",replacement = "",x = str_extract(string = current_gene_nm,pattern = "\\.Y.{6}\\.[A-Z]|\\.Y.{6}"),fixed = T)
  current_gene_nm_start_X <- v_gene_nm_start_X[j]#gsub(pattern = ".",replacement = "",x = str_extract(string = current_gene_nm,pattern = "X.*?\\."),fixed = T)
  
  #if there is a match between any of the 2 ways of naming gene, make the naming uniform in both dataset
  if (!is.na(current_gene_nm_start_Y)){
    if (current_gene_nm_start_Y %in% df_Scer_expression$matchable_gene_name){
      colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_Y 
      next()
    }
  }else{
    if (current_gene_nm_start_X %in% df_Scer_expression$matchable_gene_name){
      colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_X 
      next()
    }
  }
  
    #remove "." character in the current colnames(mtx_Scer_pres_abs)[j]
  colnames(mtx_Scer_pres_abs)[j] <- ifelse(test = substr(x = colnames(mtx_Scer_pres_abs)[j],start = 1,stop = nchar(colnames(mtx_Scer_pres_abs)[j]))==".",yes = substr(x = colnames(mtx_Scer_pres_abs)[j],start = 1,stop = nchar(colnames(mtx_Scer_pres_abs)[j])-1),no = colnames(mtx_Scer_pres_abs)[j])
  
  #if there is a match between any of the 2 ways of naming gene, make the naming uniform in both dataset
    #using "X..." naming match in other columns
  if (colnames(mtx_Scer_pres_abs)[j] %in% df_Scer_expression$ORF){
    df_Scer_expression$matchable_gene_name[df_Scer_expression$ORF==(colnames(mtx_Scer_pres_abs)[j])] <- colnames(mtx_Scer_pres_abs)[j]
    next()
  }
  
  if (current_gene_nm_start_X%in%df_Scer_expression$gene_name){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_X 
    df_Scer_expression$matchable_gene_name[df_Scer_expression$gene_name==current_gene_nm_start_X] <- current_gene_nm_start_X
    next()
  }
  
  if (current_gene_nm_start_X%in%df_Scer_expression$systematic_name){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_X 
    df_Scer_expression$matchable_gene_name[df_Scer_expression$systematic_name==current_gene_nm_start_X] <- current_gene_nm_start_X
    next()
  }
  
  cond_uniq_match_X <- ((any(grepl(pattern = paste0(current_gene_nm_start_X,"."),x = unique(df_Scer_expression$Annotation_Name),fixed = T))) || (any(grepl(pattern = paste0(".",current_gene_nm_start_X),x = unique(df_Scer_expression$Annotation_Name),fixed = T))) )
  if (cond_uniq_match_X & (sum(cond_uniq_match_X)==1)){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_X 
    df_Scer_expression$matchable_gene_name[((any(grepl(pattern = paste0(current_gene_nm_start_X,"."),x = df_Scer_expression$Annotation_Name,fixed = T))) || (any(grepl(pattern = paste0(".",current_gene_nm_start_X),x = df_Scer_expression$Annotation_Name,fixed = T))) )] <- current_gene_nm_start_X
    next()
  }
  
    #using "Y..." naming match in other columns
  if (colnames(mtx_Scer_pres_abs)[j] %in% df_Scer_expression$ORF){
    df_Scer_expression$matchable_gene_name[df_Scer_expression$ORF==(colnames(mtx_Scer_pres_abs)[j])] <- colnames(mtx_Scer_pres_abs)[j]
    next()
  }
  
  if (current_gene_nm_start_Y%in%df_Scer_expression$gene_name){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_Y 
    df_Scer_expression$matchable_gene_name[df_Scer_expression$gene_name==current_gene_nm_start_Y] <- current_gene_nm_start_Y
    next()
  }
  
  if (current_gene_nm_start_Y%in%df_Scer_expression$systematic_name){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_Y 
    df_Scer_expression$matchable_gene_name[df_Scer_expression$systematic_name==current_gene_nm_start_Y] <- current_gene_nm_start_Y
    next()
  }
  
  cond_uniq_match_Y <- ((any(grepl(pattern = paste0(current_gene_nm_start_Y,"."),x = unique(df_Scer_expression$Annotation_Name),fixed = T))) || (any(grepl(pattern = paste0(".",current_gene_nm_start_Y),x = unique(df_Scer_expression$Annotation_Name),fixed = T))) )
  if (cond_uniq_match_Y & (sum(cond_uniq_match_Y)==1)){
    colnames(mtx_Scer_pres_abs)[j] <- current_gene_nm_start_Y 
    df_Scer_expression$matchable_gene_name[((any(grepl(pattern = paste0(current_gene_nm_start_Y,"."),x = df_Scer_expression$Annotation_Name,fixed = T))) || (any(grepl(pattern = paste0(".",current_gene_nm_start_Y),x = df_Scer_expression$Annotation_Name,fixed = T))) )] <- current_gene_nm_start_Y
    next()
  }
}

  #Create the matrix of expression with genes for with presence/absence data
mtx_Scer_expression_for_genes_with_pres_abs_data <- reshape2::acast(subset(df_Scer_expression,matchable_gene_name%in%colnames(mtx_Scer_pres_abs)), Strain ~ matchable_gene_name, value.var = "tpm",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)
mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data <- t(mtx_Scer_expression_for_genes_with_pres_abs_data)
#remove non-expressed genes
dim(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)
mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data <- mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data[rowSums(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)!=0,]
dim(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)

proportion_core_genes_in_Scer_pangenome <- sum(rowSums(t(mtx_Scer_pres_abs),na.rm = T)/ncol(t(mtx_Scer_pres_abs))==1)/nrow(t(mtx_Scer_pres_abs))
proportion_core_genes_in_Scer_pangenome
proportion_accessory_genes_in_Scer_pangenome <- 1 - proportion_core_genes_in_Scer_pangenome
proportion_accessory_genes_in_Scer_pangenome
mtx_Scer_pres_abs_filtered_for_genes_with_expr_data <- mtx_Scer_pres_abs[,colnames(mtx_Scer_pres_abs)[colnames(mtx_Scer_pres_abs)%in%df_Scer_expression$matchable_gene_name]]
mtx_TRANSPOSE_Scer_pres_abs_filtered_for_genes_with_expr_data <- t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data)
#make sure to keep the same genes as in mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data
dim(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)
mtx_TRANSPOSE_Scer_pres_abs_filtered_for_genes_with_expr_data <- mtx_TRANSPOSE_Scer_pres_abs_filtered_for_genes_with_expr_data[rownames(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data),]
dim(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)

  #pca TRANSPOSED Sc expression
#Perform the PCA for both conditions
Scer_TRANSPOSE_gene_expr_pca <- prcomp(mtx_TRANSPOSE_Scer_expression_for_genes_with_pres_abs_data)
#Visualize PCA
Scer_TRANSPOSE_gene_expr_pc_eigenvalues <- Scer_TRANSPOSE_gene_expr_pca$sdev^2

df_Scer_TRANSPOSE_gene_expr_pca <- as.data.frame(as_tibble(Scer_TRANSPOSE_gene_expr_pca$x, rownames = "species"))
df_Scer_TRANSPOSE_gene_expr_pca$species_initials <- vapply(X = df_Scer_TRANSPOSE_gene_expr_pca$species,FUN = get_initials,FUN.VALUE = "")
df_Scer_TRANSPOSE_gene_expr_pca$mixture <- v_species_smpl_lbl_to_mixture[df_Scer_TRANSPOSE_gene_expr_pca$species]
#Eigenvalues
Scer_TRANSPOSE_gene_expr_pc_eigenvalues <- Scer_TRANSPOSE_gene_expr_pca$sdev^2
df_Scer_TRANSPOSE_gene_expr_pc_eigenvalues <- tibble(PC = factor(1:length(Scer_TRANSPOSE_gene_expr_pc_eigenvalues)), 
                                                            variance = Scer_TRANSPOSE_gene_expr_pc_eigenvalues) %>% 
  #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
  #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

nb_pcs_recapitulating_Scer_TRANSPOSE_gene_expr_pca <- as.integer(subset(df_Scer_TRANSPOSE_gene_expr_pc_eigenvalues,explained_variance_cum>=0.99)$PC)[1]
#Visualize species aggregate gene expression PCA eigen values
df_Scer_TRANSPOSE_gene_expr_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = explained_variance)) +
  geom_line(aes(y = explained_variance_cum, group = 1)) + 
  geom_point(aes(y = explained_variance_cum)) +
  geom_hline(yintercept = 0.99,col="red",lty=2)+
  geom_vline(xintercept = nb_pcs_recapitulating_Scer_TRANSPOSE_gene_expr_pca,col="red",lty=2)+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  labs(x = "Principal component", y = "Fraction variance explained")
ggsave(filename = "PCA_Scer_TRANSPOSE_gene_expr_Cumul_Explained_variance.png", path=output_workspace, width = 25/2.54, height = 15/2.54)

#Visualize species aggregate gene expression PCA
df_Scer_TRANSPOSE_gene_expr_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(size = 2, stroke = 2) + 
  ggrepel::geom_text_repel(aes(label = mixture)) #geom_text
ggsave(filename = "PCA_Scer_TRANSPOSE_gene_expr.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

  #PCs multiple linear regressions with presence/absence
    #Based on multiple linear regression (MLR)
v_Scer_gene_expr_Transpose_pcs_var_ratio <- df_Scer_TRANSPOSE_gene_expr_pc_eigenvalues$explained_variance
v_Scer_gene_expr_Transpose_pc_vs_Pres_abs_R2 <- rep(NA,length(v_Scer_gene_expr_Transpose_pcs_var_ratio))
for (the_pc_id in 1:length(v_Scer_gene_expr_Transpose_pcs_var_ratio)){
  v_Scer_gene_expr_Transpose_pc_vs_Pres_abs_R2[the_pc_id] <- summary(lm(Scer_TRANSPOSE_gene_expr_pca$x[,the_pc_id]~mtx_TRANSPOSE_Scer_pres_abs_filtered_for_genes_with_expr_data))$r.squared
}
var_ratio_of_Scer_expression_explained_by_presence_absence_MLR <- sum(v_Scer_gene_expr_Transpose_pcs_var_ratio*v_Scer_gene_expr_Transpose_pc_vs_Pres_abs_R2)
var_ratio_of_Scer_expression_explained_by_presence_absence_MLR
      #bootstrap var ratio of S. cerevisiae expression_explained_by_presence_absence
nb_subsamples <- 999
v_BOOTSTRAP_var_ratio_of_Scer_expression_explained_by_presence_absence_MLR <- NULL
for (the_current_permutation_id in 1:nb_subsamples){
  v_current_permutation_gene_indices <- sample(x = 1:nrow(Scer_TRANSPOSE_gene_expr_pca$x),size = 100,replace = F)
  v_BOOTSTRAP_S_cer_gene_expr_Transpose_pc_vs_Pres_abs_R2 <- rep(NA,length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB))
  for (the_pc_id in 1:length(v_species_aggregate_gene_expr_Transpose_pcs_var_ratio_YNB)){
    v_BOOTSTRAP_S_cer_gene_expr_Transpose_pc_vs_Pres_abs_R2[the_pc_id] <- summary(lm(Scer_TRANSPOSE_gene_expr_pca$x[v_current_permutation_gene_indices,the_pc_id]~mtx_TRANSPOSE_Scer_pres_abs_filtered_for_genes_with_expr_data[v_current_permutation_gene_indices,]))$r.squared
  }
  v_BOOTSTRAP_var_ratio_of_Scer_expression_explained_by_presence_absence_MLR <- c(v_BOOTSTRAP_var_ratio_of_Scer_expression_explained_by_presence_absence_MLR,sum(v_Scer_gene_expr_Transpose_pcs_var_ratio*v_BOOTSTRAP_S_cer_gene_expr_Transpose_pc_vs_Pres_abs_R2))
  if (the_current_permutation_id%%5==0){
    print(the_current_permutation_id)
    print(v_BOOTSTRAP_S_cer_gene_expr_Transpose_pc_vs_Pres_abs_R2)
  }
}
summary(v_BOOTSTRAP_var_ratio_of_Scer_expression_explained_by_presence_absence_MLR);sd(v_BOOTSTRAP_var_ratio_of_Scer_expression_explained_by_presence_absence_MLR)

#S. cerevisiae PCA (REGULAR FORMAT NOT TRANSPOSED)
mtx_Scer_expression_all_genes <- reshape2::acast(df_Scer_expression, Strain ~ matchable_gene_name, value.var = "tpm",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)
  #Perform the PCA for both conditions
Scer_regular_gene_expr_pca <- prcomp(mtx_Scer_expression_all_genes)
  #Visualize PCA
Scer_regular_gene_expr_pc_eigenvalues <- Scer_regular_gene_expr_pca$sdev^2

df_Scer_regular_gene_expr_pca <- as.data.frame(as_tibble(Scer_regular_gene_expr_pca$x, rownames = "species"))
df_Scer_regular_gene_expr_pca$species_initials <- vapply(X = df_Scer_regular_gene_expr_pca$species,FUN = get_initials,FUN.VALUE = "")
df_Scer_regular_gene_expr_pca$mixture <- v_species_smpl_lbl_to_mixture[df_Scer_regular_gene_expr_pca$species]
  #Eigenvalues
Scer_regular_gene_expr_pc_eigenvalues <- Scer_regular_gene_expr_pca$sdev^2
df_Scer_regular_gene_expr_pc_eigenvalues <- tibble(PC = factor(1:length(Scer_regular_gene_expr_pc_eigenvalues)), 
                                                            variance = Scer_regular_gene_expr_pc_eigenvalues) %>% 
  #add a new column with the percent variance
  mutate(explained_variance = variance/sum(variance)) %>% 
  #add another column with the cumulative variance explained
  mutate(explained_variance_cum = cumsum(explained_variance)) %>%
  as.data.frame()

nb_pcs_recapitulating_Scer_regular_gene_expr_pca <- as.integer(subset(df_Scer_regular_gene_expr_pc_eigenvalues,explained_variance_cum>=0.99)$PC)[1]

#save Scer strains regular gene expression matrix PCA
write.table(x=df_Scer_regular_gene_expr_pc_eigenvalues,file = paste0(output_workspace,"Table_df_Scer_regular_gene_expr_pc_eigenvalues.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=Scer_regular_gene_expr_pca$x,file = paste0(output_workspace,"mtx_Scer_regular_gene_expr_pca.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=Scer_regular_gene_expr_pca$x[,1:nb_pcs_recapitulating_Scer_regular_gene_expr_pca],file = paste0(output_workspace,"mtx_Scer_regular_gene_expr_pca_RETAINED_PCs.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

  #Visualize species aggregate gene expression PCA eigen values
df_Scer_regular_gene_expr_pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = explained_variance)) +
  geom_line(aes(y = explained_variance_cum, group = 1)) + 
  geom_point(aes(y = explained_variance_cum)) +
  geom_hline(yintercept = 0.99,col="red",lty=2)+
  geom_vline(xintercept = nb_pcs_recapitulating_Scer_regular_gene_expr_pca,col="red",lty=2)+
  scale_y_continuous(breaks = seq(0,1,0.2),limits = c(0,1))+
  labs(x = "Principal component", y = "Fraction variance explained")
ggsave(filename = "PCA_Scer_regular_gene_expr_Cumul_Explained_variance.png", path=output_workspace, width = 25/2.54, height = 15/2.54)

  #Visualize species aggregate gene expression PCA
df_Scer_regular_gene_expr_pca %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(size = 2, stroke = 2) + 
  ggrepel::geom_text_repel(aes(label = mixture)) #geom_text
ggsave(filename = "PCA_Scer_regular_gene_expr.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

  #Visualize top 10 ORTHOGROUPS (or "gene") with highest loadings in species aggregate gene expression PCA's pc1 and pc2
Scer_regular_gene_expr_pc_loadings <- Scer_regular_gene_expr_pca$rotation %>% 
  as_tibble(rownames = "gene") %>%
  as.data.frame()

v_Scer_regular_gene_expr_pc_eigenvalues <- Scer_regular_gene_expr_pc_eigenvalues
names(v_Scer_regular_gene_expr_pc_eigenvalues) <- paste0("PC",1:length(Scer_regular_gene_expr_pc_eigenvalues))

  #as pc1/pc2 loading dataframe
df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings <- Scer_regular_gene_expr_pc_loadings %>% 
  #select only the PCs we are interested in
  dplyr::select(gene, PC1, PC2) %>%
  #convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  #for each PC
  dplyr::group_by(PC) %>%
  as.data.frame()

df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$the_princ_comp_explained_var <- v_Scer_regular_gene_expr_pc_eigenvalues[df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$PC]
df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$normalized_loading <- df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$loading*df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$the_princ_comp_explained_var

df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings <- df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings %>% 
  #arrange by descending order of loading
  dplyr::arrange(desc(abs(normalized_loading))) 

Scer_regular_gene_expr_top10_gene_orthogroup <- df_sorted_Scer_regular_gene_expr_pc1_and_pc2_loadings$gene[1:10]

Scer_regular_gene_expr_top10_gene_orthogroup 
  #Functions Scer_regular_gene_expr_top10_gene_orthogroup:
Scer_regular_gene_expr_top10_gene_orthogroup_functions <- c("YMR251W-A stress response protein (hyperosmolarity; inositol-deficiency-induced Calcium sensitivity)","YPR036W-A sporulation","YPR036W-A sporulation","YFL014W Plasma membrane protein involved in maintaining membrane organization; involved in maintaining organization during stress conditions; induced by heat shock, oxidative stress, osmostress, stationary phase, glucose depletion, oleate and alcohol","YFL014W Plasma membrane protein involved in maintaining membrane organization; involved in maintaining organization during stress conditions; induced by heat shock, oxidative stress, osmostress, stationary phase, glucose depletion, oleate and alcohol","YMR251W-A stress response protein (hyperosmolarity; inositol-deficiency-induced Calcium sensitivity)","YGR192C TDH3 Triose-phosphate DeHydrogenase (glycolysis and gluconeogenesis)","YLR327C Translation Machinery Associated protein involved in stress response","YLR167W Ribosomal Protein of the Small subunit","YGR192C TDH3 Triose-phosphate DeHydrogenase (glycolysis and gluconeogenesis)")
Scer_regular_gene_expr_top10_gene_orthogroup_functions

  #visualize sample PC1/PC2 top 10 orthogroup pc loading
Scer_regular_gene_expr_top_loadings <- Scer_regular_gene_expr_pc_loadings %>% 
  dplyr::filter(gene %in% Scer_regular_gene_expr_top10_gene_orthogroup)

ggplot(data = Scer_regular_gene_expr_top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               col = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = gene),
                           nudge_y = 0.005, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))+
  labs(x = "PC1", y = "PC2")
ggsave(filename = "Top10_pc1_pc2_normalized_loadings_PCA_Scer_regular_gene_expr.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

#**********************************************************************Fitness data and distances*************************************************************************************************#

#S. cerevisiae strains Fitness_distance ~ Expression_distance + Gene_content_distance + Phylogenetic_distance
  #Import S.cer. fitness data
df_Scer_strains_fitness <- read.csv2(file = paste0(output_workspace,"Sc_pheno_35Conditions_NormalizedByYPD.txt"),sep = "\t",header = T,stringsAsFactors = FALSE)
df_Scer_strains_fitness$YPD40 <- as.numeric(df_Scer_strains_fitness$YPD40)
  #Get the S.cer strains fitness distance
mtx_Scer_strains_fitness <- as.matrix(df_Scer_strains_fitness$YPD40)
rownames(mtx_Scer_strains_fitness) <- df_Scer_strains_fitness$X
mtx_dist_Scer_strains_fitness <- as.matrix(dist(x = mtx_Scer_strains_fitness,method = "euclidean"))
#remove(mtx_Scer_strains_fitness)
dim(mtx_dist_Scer_strains_fitness)
  #Get the expression, gene_content and phylo distances
mtx_dist_expression_Scer_strains <- as.matrix(1/as.matrix(cor(t(mtx_Scer_expression_all_genes))^2))#as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD,method = "euclidean"))
mtx_dist_gene_content_Scer_strains <- as.matrix(vegan::vegdist(mtx_Scer_pres_abs,method = "gower",na.rm = T))
mtx_phylo_dist_Scer_strains <- as.matrix(read.csv2(file = paste0(output_workspace,"Scer_1011_strains_distanceMatrixBasedOnSNPS.tab"),sep = "\t",header = T,row.names = 1,stringsAsFactors = FALSE))
mtx_phylo_dist_Scer_strains <- apply(mtx_phylo_dist_Scer_strains, 2, as.numeric)
rownames(mtx_phylo_dist_Scer_strains) <- colnames(mtx_phylo_dist_Scer_strains)
#function to get the upper triangle part of a matrix without the diagonal and in a vector format
get_v_uppr_tri_no_diag_from_mtx <- function(the_mtx){
  the_mtx[upper.tri(the_mtx, diag = FALSE)]
}

  #Find the strains that are present in all 4 matrices
v_Scer_matrices_strain_intersect <- intersect(rownames(mtx_phylo_dist_Scer_strains),intersect(rownames(mtx_dist_gene_content_Scer_strains),intersect(rownames(mtx_dist_expression_Scer_strains),rownames(mtx_dist_Scer_strains_fitness))))

#Find the expression distances while controlling for gene content (only on genes shared by species pairs)
mtx_gene_content_with_NA_for_absence_in_Scer_strains <- mtx_Scer_pres_abs_filtered_for_genes_with_expr_data[v_Scer_matrices_strain_intersect,]
mtx_gene_content_with_NA_for_absence_in_Scer_strains[mtx_gene_content_with_NA_for_absence_in_Scer_strains==0] <- NA
mtx_pairwise_expression_similarity_based_on_present_genes_in_Scer_strains <- as.matrix(cor(t(mtx_Scer_expression_for_genes_with_pres_abs_data*mtx_gene_content_with_NA_for_absence_in_Scer_strains[v_Scer_matrices_strain_intersect,]),use = "pairwise.complete.obs")^2)

#apply the fct to mtx_dist_Scer_strains_fitness
v_dist_Scer_strains_fitness <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_Scer_strains_fitness[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect])
(length(v_dist_Scer_strains_fitness)==(length(v_Scer_matrices_strain_intersect)*length(v_Scer_matrices_strain_intersect)/2)-(length(v_Scer_matrices_strain_intersect)/2))

  #VarPart of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
    #create the dataframe of explanatory variable of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
df_expl_vars_Scer_distances_varpart <- data.frame(Expression_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_expression_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect]),Gene_content_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_gene_content_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect]),CoreGenome_SNPs_Phylogenetic_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_phylo_dist_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect]),stringsAsFactors = T)
      #df_expl_vars_Scer_distances_varpart to long format
#df_expl_vars_Scer_distances_varpart_long <- df_expl_vars_Scer_distances_varpart %>%
#  pivot_longer(
#    cols = c(Expression_distance,Gene_content_distance,CoreGenome_SNPs_Phylogenetic_distance), 
#    names_to = "Explanatory_variable", 
#    values_to = "value") %>% as.data.frame()
df_Scer_distances_varpart <- df_expl_vars_Scer_distances_varpart
df_Scer_distances_varpart$Fitness_distance <- v_dist_Scer_strains_fitness
remove(df_expl_vars_Scer_distances_varpart)
Scer_distances_varpart <- vegan::varpart(as.matrix(df_Scer_distances_varpart$Fitness_distance), as.matrix(df_Scer_distances_varpart$Expression_distance),as.matrix(df_Scer_distances_varpart$Gene_content_distance), as.matrix(df_Scer_distances_varpart$CoreGenome_SNPs_Phylogenetic_distance),scale = F)
Scer_distances_varpart
summary(lm(data = df_Scer_distances_varpart,formula = Fitness_distance~Expression_distance+Gene_content_distance+CoreGenome_SNPs_Phylogenetic_distance))

#Using fitness variation vs each predictor covariance matrix (covariance-based VARPART method)
  #Expression
#summary(lm(formula = "fitness~.",data = cbind(data.frame(fitness=df_Scer_strains_fitness$YPD40[df_Scer_strains_fitness$X%in%v_Scer_matrices_strain_intersect]),as.data.frame(cov(t(mtx_Scer_expression_for_genes_with_pres_abs_data[v_Scer_matrices_strain_intersect,apply(X = mtx_Scer_expression_for_genes_with_pres_abs_data,MARGIN = 2,FUN = sd)!=0]))/sum(apply(X = mtx_Scer_expression_for_genes_with_pres_abs_data,MARGIN = 2,FUN = sd)!=0) ) )))
  #Phylogenetic (Core genome SNPs and divergence)
#summary(lm(formula = "fitness~.",data = cbind(data.frame(fitness=df_Scer_strains_fitness$YPD40[df_Scer_strains_fitness$X%in%v_Scer_matrices_strain_intersect]),as.data.frame(cov(t(mtx_phylo_dist_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect]))/length(v_Scer_matrices_strain_intersect) ) ) ))
  #Gene content 
#summary(lm(formula = "fitness~.",data = cbind(data.frame(fitness=df_Scer_strains_fitness$YPD40[df_Scer_strains_fitness$X%in%v_Scer_matrices_strain_intersect]),as.data.frame(cov(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data[v_Scer_matrices_strain_intersect,]),use = "pairwise.complete.obs")/length(v_Scer_matrices_strain_intersect) ) ) ))
#SEE REML for results

#Yeast species Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
  #Import yeast species fitness data
df_yeast_species_fitness <- read.csv2(file = paste0(output_workspace,"y1000p_Growth_data.csv"),sep = ",",header = T,stringsAsFactors = FALSE)
df_yeast_species_fitness$Species <- ifelse(test = substr(x = df_yeast_species_fitness$Species,start = nchar(df_yeast_species_fitness$Species),stop = nchar(df_yeast_species_fitness$Species) )==" ",yes = substr(x = df_yeast_species_fitness$Species,start = 1,stop = nchar(df_yeast_species_fitness$Species)-1),no = df_yeast_species_fitness$Species )
df_yeast_species_fitness$YPD <- as.numeric(df_yeast_species_fitness$Glucose)
df_yeast_species_fitness$YNB <- as.numeric(df_yeast_species_fitness$No.Carbon)
  #Get the S.cer strains fitness distance
    #YPD
mtx_yeast_species_fitness_YPD <- as.matrix(df_yeast_species_fitness$YPD)
rownames(mtx_yeast_species_fitness_YPD) <- df_yeast_species_fitness$Species
mtx_dist_yeast_species_fitness_YPD <- as.matrix(dist(x = mtx_yeast_species_fitness_YPD,method = "euclidean"))
#remove(mtx_yeast_species_fitness_YPD)
dim(mtx_dist_yeast_species_fitness_YPD)
  #make sure the species are named the same way in the phylogenetic matrix
mtx_phylo_dist_the_sps_YPD <- mtx_phylo_dist_the_21_sps_YPD
rownames(mtx_phylo_dist_the_sps_YPD) <- gsub(pattern = "_",replacement = " ",x = rownames(mtx_phylo_dist_the_sps_YPD),fixed = T)
#apply the get_v_uppr_tri_no_diag_from_mtx to mtx_dist_yeast_species_fitness_YPD
v_dist_fitness_21_spcs_YPD <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_yeast_species_fitness_YPD[v_lst_species_smpl_lbl_YPD,v_lst_species_smpl_lbl_YPD])
(length(v_dist_fitness_21_spcs_YPD)==(length(v_lst_species_smpl_lbl_YPD)*length(v_lst_species_smpl_lbl_YPD)/2)-(length(v_lst_species_smpl_lbl_YPD)/2))

#VarPart of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
#create the dataframe of explanatory variable of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
#df_expl_vars_Yeast_21_spcs_distances_varpart_YPD <- data.frame(Expression_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_aggregate_expr_YPD[v_lst_species_smpl_lbl_YPD,v_lst_species_smpl_lbl_YPD]),Gene_content_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_pres_abs_YPD[v_lst_species_smpl_lbl_YPD,v_lst_species_smpl_lbl_YPD]),CoreGenome_SNPs_Phylogenetic_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_phylo_dist_the_21_sps_YPD[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD),gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD)]),stringsAsFactors = T)
#df_Yeast_21_spcs_distances_varpart_YPD <- df_expl_vars_Yeast_21_spcs_distances_varpart_YPD
#df_Yeast_21_spcs_distances_varpart_YPD$Fitness_distance <- v_dist_fitness_21_spcs_YPD
#remove(df_expl_vars_Yeast_21_spcs_distances_varpart_YPD)
#Yeast_21_spcs_distances_varpart_YPD <- vegan::varpart(as.matrix(df_Yeast_21_spcs_distances_varpart_YPD$Fitness_distance), as.matrix(df_Yeast_21_spcs_distances_varpart_YPD$Expression_distance),as.matrix(df_Yeast_21_spcs_distances_varpart_YPD$Gene_content_distance), as.matrix(df_Yeast_21_spcs_distances_varpart_YPD$CoreGenome_SNPs_Phylogenetic_distance),scale = F)
#Yeast_21_spcs_distances_varpart_YPD
#summary(lm(data = df_Yeast_21_spcs_distances_varpart_YPD,formula = Fitness_distance~Expression_distance+Gene_content_distance+CoreGenome_SNPs_Phylogenetic_distance))

#    #YNB (Fitness data are 0s)
#mtx_yeast_species_fitness_YNB <- as.matrix(df_yeast_species_fitness$YNB)
#rownames(mtx_yeast_species_fitness_YNB) <- df_yeast_species_fitness$Species
#mtx_dist_yeast_species_fitness_YNB <- as.matrix(dist(x = mtx_yeast_species_fitness_YNB,method = "euclidean"))
#remove(mtx_yeast_species_fitness_YNB)
#dim(mtx_dist_yeast_species_fitness_YNB)
##make sure the species are named the same way in the phylogenetic matrix
#mtx_phylo_dist_the_sps_YNB <- mtx_phylo_dist_the_21_sps_YNB
#rownames(mtx_phylo_dist_the_sps_YNB) <- gsub(pattern = "_",replacement = " ",x = rownames(mtx_phylo_dist_the_sps_YNB),fixed = T)
##apply the get_v_uppr_tri_no_diag_from_mtx to mtx_dist_yeast_species_fitness_YNB
#v_dist_fitness_21_spcs_YNB <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_yeast_species_fitness_YNB[v_lst_species_smpl_lbl_YNB,v_lst_species_smpl_lbl_YNB])
#(length(v_dist_fitness_21_spcs_YNB)==(length(v_lst_species_smpl_lbl_YNB)*length(v_lst_species_smpl_lbl_YNB)/2)-(length(v_lst_species_smpl_lbl_YNB)/2))
#
##VarPart of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
##create the dataframe of explanatory variable of Fitness_distance ~ Expression_distance + Gene_content_distance + CoreGenome_SNPs_Phylogenetic_distance
#df_expl_vars_Yeast_21_spcs_distances_varpart_YNB <- data.frame(Expression_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_aggregate_expr_YNB[v_lst_species_smpl_lbl_YNB,v_lst_species_smpl_lbl_YNB]),Gene_content_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_pres_abs_YNB[v_lst_species_smpl_lbl_YNB,v_lst_species_smpl_lbl_YNB]),CoreGenome_SNPs_Phylogenetic_distance=get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_phylo_dist_the_21_sps_YNB[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YNB),gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YNB)]),stringsAsFactors = T)
#df_Yeast_21_spcs_distances_varpart_YNB <- df_expl_vars_Yeast_21_spcs_distances_varpart_YNB
#df_Yeast_21_spcs_distances_varpart_YNB$Fitness_distance <- v_dist_fitness_21_spcs_YNB
#remove(df_expl_vars_Yeast_21_spcs_distances_varpart_YNB)
#Yeast_21_spcs_distances_varpart_YNB <- vegan::varpart(as.matrix(df_Yeast_21_spcs_distances_varpart_YNB$Fitness_distance), as.matrix(df_Yeast_21_spcs_distances_varpart_YNB$Expression_distance),as.matrix(df_Yeast_21_spcs_distances_varpart_YNB$Gene_content_distance), as.matrix(df_Yeast_21_spcs_distances_varpart_YNB$CoreGenome_SNPs_Phylogenetic_distance),scale = F)
#Yeast_21_spcs_distances_varpart_YNB
#summary(lm(data = df_Yeast_21_spcs_distances_varpart_YNB,formula = Fitness_distance~Expression_distance+Gene_content_distance+CoreGenome_SNPs_Phylogenetic_distance))


#Expression plasticity YPD --> YNB
  #Based on distances (which include the impact of differentially expressed genes)
summary(lm(as.vector(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[v_lst_species_smpl_lbl_YNB, intersect(colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD),colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB))])~as.vector(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[v_lst_species_smpl_lbl_YNB, intersect(colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD),colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB))])))
  #Based on shared expressed genes transcript level
summary(lm(get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_aggregate_expr_YPD[v_lst_species_smpl_lbl_YNB,v_lst_species_smpl_lbl_YNB])~get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_aggregate_expr_YNB[v_lst_species_smpl_lbl_YNB,v_lst_species_smpl_lbl_YNB])))

#Save all matrices required for GREML Varpart (MAKE SURE THE ROW ORDER IS THE SAME!!!)
  #Scer strains
    #Fitness mtx
write.table(x=unname(mtx_Scer_strains_fitness[v_Scer_matrices_strain_intersect,]),file = paste0(output_workspace,"lst_Scer_strains_Fitness.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Expression mtx
write.table(x=mtx_Scer_expression_all_genes[v_Scer_matrices_strain_intersect,],file = paste0(output_workspace,"mtx_Scer_strains_Expression.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Gene presense/absence
write.table(x=mtx_Scer_pres_abs[v_Scer_matrices_strain_intersect,],file = paste0(output_workspace,"mtx_Scer_strains_Gene_Presence_Absence.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Phylo distance kernel
write.table(x=mtx_phylo_dist_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect],file = paste0(output_workspace,"mtx_Scer_strains_Phylo_dist_Kernel.csv"),sep = ",",na = "NA",row.names = F,col.names = F)

  #Yeast 21 spcs YPD (YNB fitness data have no variance because fitness = 0 for all species in the Science paper)
    #Fitness mtx
write.table(x=unname(mtx_yeast_species_fitness_YPD[v_lst_species_smpl_lbl_YPD,]),file = paste0(output_workspace,"lst_Yeast_21_spcs_Fitness_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Expression mtx
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[v_lst_species_smpl_lbl_YPD,],file = paste0(output_workspace,"mtx_Yeast_21_spcs_Expression_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Gene presense/absence
write.table(x=mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[v_lst_species_smpl_lbl_YPD,],file = paste0(output_workspace,"mtx_Yeast_21_spcs_Gene_Presence_Absence_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Phylo distance kernel
write.table(x=mtx_phylo_dist_the_21_sps_YPD[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD),gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD)],file = paste0(output_workspace,"mtx_Yeast_21_spcs_Phylo_dist_Kernel_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)

    #vcv matrix
newick_string <- "(Schizosaccharomyces_pombe:1.32455,(Lipomyces_kononenkoae:0.75795,((Blastobotrys_peoriensis:0.71357,Sugiyamaella_smithiae:0.67961):0.27217[100],((Saccharomycopsis_amapae:0.87613,((Wickerhamomyces_anomalus:0.26769,(Candida_ulmi:0.44249,(Wickerhamomyces_strasburgensis:0.28446,(Cyberlindnera_mrakii:0.33120,Cyberlindnera_japonica:0.33370):0.06143[100]):0.07477[100]):0.06937[100]):0.29750[100],((Hanseniaspora_opuntiae:0.13825,Hanseniaspora_nectarophila:0.16210):1.57101[100],Saccharomyces_cerevisiae:0.72534):0.37896[100]):0.07611[100]):0.03832[100],(((Nakazawaea_peltata:0.42461,Nakazawaea_wyomingensis:0.31925):0.40512[100],(Ambrosiozyma_oregonensis:0.48889,Saturnispora_besseyi:0.94816):0.49401[100]):0.10894[100],((Meyerozyma_caribbica:0.55152,Candida_blattariae:0.74028):0.01876[100],(Candida_dosseyi:0.69857,((Schwanniomyces_capriottii:0.06750,Schwanniomyces_polymorphus_var._polymorphus:0.15034):0.32983[100],Kodamaea_ohmeri:0.76246):0.04910[88]):0.02024[100]):0.44007[100]):0.04475[100]):0.38294[100]):0.18123[100]):0.07902):0.00000;"
tree <- read.tree(text = newick_string)
vcv_matrix <- vcv(tree)
vcv_matrix <- vcv_matrix[gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD),gsub(pattern = " ",replacement = "_",x = v_lst_species_smpl_lbl_YPD)]
#save vcv matric
write.table(x=vcv_matrix,file = paste0(output_workspace,"Phylo_vcv_kernel_21_spcs_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
    #Regular Expression PCA 
write.table(x=species_aggregate_gene_expr_pca_YPD$x[,1:nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YPD],file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_YPD_RETAINED_PCs.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)
write.table(x=species_aggregate_gene_expr_pca_YNB$x[,1:nb_pcs_recapitulating_species_aggregate_gene_expr_pca_YNB],file = paste0(output_workspace,"mtx_species_aggregate_gene_expr_pca_YNB_RETAINED_PCs.tsv"),sep = "\t",na = "NA",row.names = T,col.names = T)

#**********************************************************************Fitness variance partitioning**************************************************************************#
#Full: F~K_E+K_pa+Phylo_dist
#P_E_pa: F~K_E+K_pa
#P_E_phylo: F~K_E+Phylo_dist
#P_pa_phylo: F~K_pa+Phylo_dist

#Fitness variance partitioning
v_comps <- c("Gene content","Expression","Phylogenetic distance")
lst_comps <- list("comp_a"= c(v_comps[1],v_comps[1]),
                  "comp_b"=c(v_comps[1],v_comps[2]),
                  "comp_c"=c(v_comps[2],v_comps[2]),
                  "comp_d"=c(v_comps[1],v_comps[3]),
                  "comp_e"=c(v_comps[1],v_comps[3]),
                  "comp_f"=c(v_comps[2],v_comps[3]),
                  "comp_g"=c(v_comps[3],v_comps[3]),
                  "residual_comp"=c("Residual","Residual"))
  #Across S cerevisiae strains (short timescale of evolution)
    #import the variance partitioning bootstrap estimates 
df_Scer_fitness_varpart_BOOTSTRAP_estimates <- read.delim(paste0(output_workspace,"Table_phenotype_variance_explained_in_the_dataset_Bootstrap_Scerevisiae_1000p_strains.tsv"))
    #filter out the bootstrap without REML convergence
df_Scer_fitness_varpart_BOOTSTRAP_estimates <- subset(df_Scer_fitness_varpart_BOOTSTRAP_estimates,(comp_Gene_content<=(1-residual_comp))&(comp_Expression<=(1-residual_comp))&(comp_Phylo_dist_ypd<=(1-residual_comp))&(residual_comp>=0))
    #summarize bootstrap results
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- df_Scer_fitness_varpart_BOOTSTRAP_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_fitness_var_explained") %>%
  group_by(Component) %>% 
  summarise(mean_variance_explained = mean(proportion_of_fitness_var_explained), median_variance_explained = median(proportion_of_fitness_var_explained), se_variance_explained = sd(proportion_of_fitness_var_explained), IQR_variance_explained = IQR(proportion_of_fitness_var_explained))
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- vapply(X = df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- ifelse(test = df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="NA",yes = NA,no=df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start)
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- vapply(X = df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- ifelse(test = df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end=="NA",yes = NA,no=df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end)
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp <- vapply(X = 1:nrow(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) paste0(unique(c(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[the_i],df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp%in%c("NA","Residual")] <- NA
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination <- df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination[df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="Residual"] <- "Residual"
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$y_min <- vapply(X = 1:nrow(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) max(0,df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]-(1*df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$y_max <- vapply(X = 1:nrow(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) min(1,df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]+(1*df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2 <- NaN
df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2[df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- v_comps[2]
        #dataframe for bottom part of the upset plot
df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats#subset(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))
for (i in 1:nrow(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats)){
  if ((!is.na(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]))&(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]!=df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[i])){
    v_current_comps <- strsplit(x = df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination[i],split = "&",fixed = T)[[1]]
    if ((length(v_current_comps)<=1)||(length(v_current_comps)>3)){
      next()
    }else if (length(v_current_comps)==2){
      #add end -> start
      current_subset1 <- df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset1$edge_end <- unname(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1)
    }else if (length(v_current_comps)==3){
      #add end -> start and middle -> end
      current_subset1 <- df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end2)
      current_subset2 <- df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset2$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset2$edge_end <- unname(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(rbind(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1),current_subset2)
    }
  }else{
    next()
  }
}

    #Upset plot Bootstrap
      #upper part
upper_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_col(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained),fill="grey60") +
  geom_pointrange(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained,ymin = y_min,ymax = y_max ))+
  ylab("Proportion of fitness variance explained") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=16),axis.text.y = element_text(size=11))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
      #bottom part
bottom_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_bottom_upset_plot_Scer_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination)) + 
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_end,levels = c(v_comps,"Residual")),group = combination))+
  geom_line(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination))+
  xlab("Components") + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_text(size=16),axis.text.x = element_text(size=9))
      #save upper and bottom parts of the upset plot
#ggsave(plot = upper_part_upsetPlot_Bootstrap_Varpart_gg, filename = "Scer_Fitness_Varpart_Bootstrap_500_13_upsetPlot.svg", path=output_workspace, width = 20/2.54, height = 15/2.54,device = svg)
#ggsave(plot = bottom_part_upsetPlot_Bootstrap_Varpart_gg,filename = "Scer_bottom_part_upsetPlot_Bootstrap_500_13_gg.svg", path=output_workspace, width = 18/2.54, height = 8/2.54,device = svg)

    #import the FULL dataset estimates
df_Scer_fitness_varpart_FULL_dataset_estimates <- read.delim(paste0(output_workspace,"Table_phenotype_variance_explained_in_the_dataset_ALL_Scerevisiae_1000p_strains.tsv"))
df_Scer_fitness_varpart_FULL_dataset_estimates_long <- df_Scer_fitness_varpart_FULL_dataset_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_fitness_var_explained")
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start <- vapply(X = df_Scer_fitness_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start <- ifelse(test = df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start=="NA",yes = NA,no=df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start)
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end <- vapply(X = df_Scer_fitness_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end <- ifelse(test = df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end=="NA",yes = NA,no=df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end)
df_Scer_fitness_varpart_FULL_dataset_estimates_long$name_comp <- vapply(X = 1:nrow(df_Scer_fitness_varpart_FULL_dataset_estimates_long),function(the_i) paste0(unique(c(df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start[the_i],df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Scer_fitness_varpart_FULL_dataset_estimates_long$name_comp[df_Scer_fitness_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Scer_fitness_varpart_FULL_dataset_estimates_long$name_comp[df_Scer_fitness_varpart_FULL_dataset_estimates_long$name_comp%in%c("NA","Residual")] <- NA
df_Scer_fitness_varpart_FULL_dataset_estimates_long$combination <- df_Scer_fitness_varpart_FULL_dataset_estimates_long$name_comp
df_Scer_fitness_varpart_FULL_dataset_estimates_long$combination[df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_start=="Residual"] <- "Residual"
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end2 <- NaN
df_Scer_fitness_varpart_FULL_dataset_estimates_long$edge_end2[df_Scer_fitness_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- v_comps[2]
##dataframe for bottom part of the upset plot
#df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats <- df_Scer_fitness_varpart_FULL_dataset_estimates_long#subset(df_Scer_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))
#for (i in 1:nrow(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats)){
#  if ((!is.na(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats$edge_start[i]))&(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats$edge_start[i]!=df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats$edge_end[i])){
#    v_current_comps <- strsplit(x = df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats$combination[i],split = "&",fixed = T)[[1]]
#    if ((length(v_current_comps)<=1)||(length(v_current_comps)>3)){
#      next()
#    }else if (length(v_current_comps)==2){
#      #add end -> start
#      current_subset1 <- df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]
#      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]$edge_end)
#      current_subset1$edge_end <- unname(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]$edge_start)
#      df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats <- rbind(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats, current_subset1)
#    }else if (length(v_current_comps)==3){
#      #add end -> start and middle -> end
#      current_subset1 <- df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]
#      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]$edge_end2)
#      current_subset2 <- df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]
#      current_subset2$edge_start <- unname(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]$edge_end)
#      current_subset2$edge_end <- unname(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats[i,]$edge_start)
#      df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats <- rbind(rbind(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats, current_subset1),current_subset2)
#    }
#  }else{
#    next()
#  }
#}
#
##Upset plot FULL dataset
##upper part
#upper_part_upsetPlot_Full_Dataset_Varpart_gg <- ggplot(data = subset(df_Scer_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))) +
#  geom_col(mapping = aes(x=reorder(combination,-proportion_of_fitness_var_explained),y=proportion_of_fitness_var_explained),fill="grey60") +
#  ylab("Proportion of fitness variance explained") + 
#  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=16),axis.text.y = element_text(size=11))+
#  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
##bottom part
#bottom_part_upsetPlot_Full_Dataset_gg <- ggplot(data = subset(df_bottom_upset_plot_Scer_fitness_varpart_FULL_estimates_summary_stats,!is.na(edge_start))) +
#  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-proportion_of_fitness_var_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination)) + 
#  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-proportion_of_fitness_var_explained),y=factor(edge_end,levels = c(v_comps,"Residual")),group = combination))+
#  geom_line(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-proportion_of_fitness_var_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination))+
#  xlab("Components") + 
#  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_text(size=16),axis.text.x = element_text(size=9))
##save upper and bottom parts of the upset plot
#ggsave(plot = upper_part_upsetPlot_Full_Dataset_Varpart_gg, filename = "Scer_Fitness_Varpart_FULL_dataset_1_912_upsetPlot.svg", path=output_workspace, width = 20/2.54, height = 15/2.54,device = svg)
#ggsave(plot = bottom_part_upsetPlot_Full_Dataset_gg,filename = "Scer_bottom_part_upsetPlot_FULL_dataset_1_912_gg.svg", path=output_workspace, width = 18/2.54, height = 8/2.54,device = svg)

    #Euler diagram
eulerr_components <- subset(df_Scer_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$proportion_of_fitness_var_explained
names(eulerr_components) <- subset(df_Scer_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$combination
#Replace any potentially negative R^2 values with 0 for plotting consistency
eulerr_components[eulerr_components < 0] <- 0

#Fit the Euler diagram to the data
fit <- eulerr::euler(eulerr_components)

svg(filename = paste0(output_workspace,"Scer_Fitness_Varpart_FULL_dataset_1_912_Euler_Diagram.svg"),width = 8.8/2.54,height = 6/2.54)
#Plot the diagram
plot(fit,
     #Set the names of the three predictors (circles)
     labels = c("Gene content", "Expression", "Phylogenetic distance", "Residual"),
     
     #Use the R^2 values as labels for the regions (rounded to 2 decimal places)
     quantities = lapply(eulerr_components, function(x) round(x, 2)),
     
     main = "Partitioning of S. cerevisiae's strains fitness variance (n=912)",
     
     #Customize colors and style
     fills = list(fill = c("#3b82f6", "#ef4444", "#22c55e"), alpha = 0.7),
     edges = list(lwd = 1.5, col = "black")
)
dev.off()

#Across Yeast species (long timescale of evolution)
#import the variance partitioning bootstrap estimates 
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates <- read.delim(paste0(output_workspace,"Table_phenotype_variance_explained_in_the_dataset_Bootstrap_Yeast_subphylum_subset.tsv"))
#filter out the bootstrap without REML convergence
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates <- subset(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates,(comp_Gene_content<=(1-residual_comp))&(comp_Expression<=(1-residual_comp))&(comp_Phylo_dist_ypd<=(1-residual_comp))&(residual_comp>=0))
#summarize bootstrap results
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_fitness_var_explained") %>%
  group_by(Component) %>% 
  summarise(mean_variance_explained = mean(proportion_of_fitness_var_explained), median_variance_explained = median(proportion_of_fitness_var_explained), se_variance_explained = sd(proportion_of_fitness_var_explained), IQR_variance_explained = IQR(proportion_of_fitness_var_explained))
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- vapply(X = df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- ifelse(test = df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="NA",yes = NA,no=df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start)
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- vapply(X = df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- ifelse(test = df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end=="NA",yes = NA,no=df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end)
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) paste0(unique(c(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[the_i],df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp%in%c("NA","Residual")] <- NA
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination <- df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$name_comp
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination[df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="Residual"] <- "Residual"
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$y_min <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) max(0,df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]-(1*df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$y_max <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) min(1,df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]+(1*df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2 <- NaN
df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2[df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- v_comps[2]
#dataframe for bottom part of the upset plot
df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats#subset(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))
for (i in 1:nrow(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats)){
  if ((!is.na(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]))&(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]!=df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[i])){
    v_current_comps <- strsplit(x = df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats$combination[i],split = "&",fixed = T)[[1]]
    if ((length(v_current_comps)<=1)||(length(v_current_comps)>3)){
      next()
    }else if (length(v_current_comps)==2){
      #add end -> start
      current_subset1 <- df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset1$edge_end <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1)
    }else if (length(v_current_comps)==3){
      #add end -> start and middle -> end
      current_subset1 <- df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end2)
      current_subset2 <- df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset2$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset2$edge_end <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(rbind(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1),current_subset2)
    }
  }else{
    next()
  }
}

#Upset plot Bootstrap
#upper part
upper_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_col(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained),fill="grey60") +
  geom_pointrange(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained,ymin = y_min,ymax = y_max ))+
  ylab("Proportion of fitness variance explained") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=16),axis.text.y = element_text(size=11))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
#bottom part
bottom_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_bottom_upset_plot_Yeast_subphylum_subset_fitness_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination)) + 
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_end,levels = c(v_comps,"Residual")),group = combination))+
  geom_line(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination))+
  xlab("Components") + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_text(size=16),axis.text.x = element_text(size=9))
#save upper and bottom parts of the upset plot
#ggsave(plot = upper_part_upsetPlot_Bootstrap_Varpart_gg, filename = "Yeast_subphylum_subset_Fitness_Varpart_Bootstrap_500_13_upsetPlot.svg", path=output_workspace, width = 20/2.54, height = 15/2.54,device = svg)
#ggsave(plot = bottom_part_upsetPlot_Bootstrap_Varpart_gg,filename = "Yeast_subphylum_subset_bottom_part_upsetPlot_Bootstrap_500_13_gg.svg", path=output_workspace, width = 18/2.54, height = 8/2.54,device = svg)

#import the FULL dataset estimates
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates <- read.delim(paste0(output_workspace,"Table_phenotype_variance_explained_in_the_dataset_ALL_Yeast_subphylum_subset.tsv"))
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long <- df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_fitness_var_explained")
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start <- vapply(X = df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start <- ifelse(test = df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start=="NA",yes = NA,no=df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start)
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end <- vapply(X = df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end <- ifelse(test = df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end=="NA",yes = NA,no=df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end)
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$name_comp <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long),function(the_i) paste0(unique(c(df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start[the_i],df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$name_comp[df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$name_comp[df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$name_comp%in%c("NA","Residual")] <- NA
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$combination <- df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$name_comp
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$combination[df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_start=="Residual"] <- "Residual"
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end2 <- NaN
df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$edge_end2[df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- v_comps[2]

#Euler diagram
eulerr_components <- subset(df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$proportion_of_fitness_var_explained
names(eulerr_components) <- subset(df_Yeast_subphylum_subset_fitness_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$combination
#Replace any potentially negative R^2 values with 0 for plotting consistency
eulerr_components[eulerr_components < 0] <- 0

#Fit the Euler diagram to the data
fit <- eulerr::euler(eulerr_components)

svg(filename = paste0(output_workspace,"Yeast_subphylum_subset_Fitness_Varpart_FULL_dataset_1_16_Euler_Diagram.svg"),width = 8.8/2.54,height = 6/2.54)
#Plot the diagram
plot(fit,
     #Set the names of the three predictors (circles)
     labels = c("Gene content", "Expression", "Phylogenetic distance", "Residual"),
     
     #Use the R^2 values as labels for the regions (rounded to 2 decimal places)
     quantities = lapply(eulerr_components, function(x) round(x, 2)),
     
     main = "Partitioning of Yeast species fitness variance (n=16)",
     
     #Customize colors and style
     fills = list(fill = c("#3b82f6", "#ef4444", "#22c55e"), alpha = 0.7),
     edges = list(lwd = 1.5, col = "black")
)
dev.off()

    #Determine the accessory orthogroups and core orthogroups BASED ON THE matrix of presence/absence
v_lst_accessory_Orthogroups_YPD <- colnames(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD)[(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD))<1)]
v_lst_core_Orthogroups_YPD <- colnames(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD)[(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD))==1)]
v_lst_accessory_Orthogroups_YNB <- colnames(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB)[(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB))<1)]
v_lst_core_Orthogroups_YNB <- colnames(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB)[(rowSums(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB))/ncol(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YNB))==1)]

    #Abundance (Mean expression) AND Dispersion (expression Mean Absolute Deviation) accessory vs core Orthogroups 21 yeast species
v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_lst_accessory_Orthogroups_YPD],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_abundance_core_Orthogroups_Yeast_21_spcs_YPD <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_lst_core_Orthogroups_YPD],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[,v_lst_accessory_Orthogroups_YNB],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_abundance_core_Orthogroups_Yeast_21_spcs_YNB <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[,v_lst_core_Orthogroups_YNB],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_lst_accessory_Orthogroups_YPD],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))
v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_lst_core_Orthogroups_YPD],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))
v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[,v_lst_accessory_Orthogroups_YNB],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))
v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB <- apply(X = mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[,v_lst_core_Orthogroups_YNB],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))

df_expr_abundance_pangenome_Yeast_21_spcs <- data.frame(Orthogroup=c(names(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),names(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),names(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB),names(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB),names(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),names(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD),names(v_abundance_core_Orthogroups_Yeast_21_spcs_YNB),names(v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB) ),
                                                        Media=rep(c("YPD","YNB","YPD","YNB"),times=c(length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YNB)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB) )),
                                                        Metric=rep(c("Abundance","Dispersion","Abundance","Dispersion","Abundance","Dispersion","Abundance","Dispersion"),times=c(length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB),length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YNB),length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB) )),
                                                        Pangenomic_component=rep(c("Accessory","Core"),times=c(length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD)+length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD)+length(v_abundance_core_Orthogroups_Yeast_21_spcs_YNB)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB) )),
                                                        Value=c(unname(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),unname(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),unname(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YNB),unname(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YNB),unname(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),unname(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD),unname(v_abundance_core_Orthogroups_Yeast_21_spcs_YNB),unname(v_dispersion_core_Orthogroups_Yeast_21_spcs_YNB) ),stringsAsFactors = F)
ggplot(data = df_expr_abundance_pangenome_Yeast_21_spcs,aes(x=as.factor(Pangenomic_component),y = Value)) + geom_violin(fill="grey") + geom_boxplot(width=0.075) + xlab("Pangenomic component") + ylab("Value") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12)) + stat_compare_means(method = "wilcox") + facet_grid(Metric~Media,scales="free_y")
ggsave(filename = "Expression_Abundance_AND_Dispersion_Accessory_vs_Core_Orthogroup_21_Yeast_spcs_YPD.png", path=output_workspace, width = 25, height = 20, units = "cm",dpi = 1200)
ggsave(filename = "Expression_Abundance_AND_Dispersion_Accessory_vs_Core_Orthogroup_21_Yeast_spcs_YPD.eps", path=output_workspace, width = 25, height = 20, units = "cm",dpi = 1200,device=cairo_ps)

  #Save expression matrices (accessory + core) with the same species order as in fitness vector
#Expression mtx
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[v_lst_species_smpl_lbl_YPD,v_lst_accessory_Orthogroups_YPD],file = paste0(output_workspace,"mtx_Yeast_21_spcs_ACCESSORY_Orthogroup_aggr_Expression_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[v_lst_species_smpl_lbl_YPD,v_lst_core_Orthogroups_YPD],file = paste0(output_workspace,"mtx_Yeast_21_spcs_CORE_Orthogroup_aggr_Expression_YPD.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[v_lst_species_smpl_lbl_YNB,v_lst_accessory_Orthogroups_YNB],file = paste0(output_workspace,"mtx_Yeast_21_spcs_ACCESSORY_Orthogroup_aggr_Expression_YNB.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=mtx_log2_1p_species_aggregate_gene_expr_FPKM_YNB[v_lst_species_smpl_lbl_YNB,v_lst_core_Orthogroups_YNB],file = paste0(output_workspace,"mtx_Yeast_21_spcs_CORE_Orthogroup_aggr_Expression_YNB.csv"),sep = ",",na = "NA",row.names = F,col.names = F)

write.table(x=mtx_Scer_expression_all_genes[,((colSums(as.matrix(mtx_Scer_expression_all_genes>0,nrow=nrow(mtx_Scer_expression_all_genes)),na.rm=T)/nrow(as.matrix(mtx_Scer_expression_all_genes>0,nrow=nrow(mtx_Scer_expression_all_genes))))<1)],file = paste0(output_workspace,"mtx_Scer_ACCESSORY_genes_Expression.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=mtx_Scer_expression_all_genes[,((colSums(as.matrix(mtx_Scer_expression_all_genes>0,nrow=nrow(mtx_Scer_expression_all_genes)),na.rm=T)/nrow(as.matrix(mtx_Scer_expression_all_genes>0,nrow=nrow(mtx_Scer_expression_all_genes))))==1)],file = paste0(output_workspace,"mtx_Scer_CORE_genes_Expression.csv"),sep = ",",na = "NA",row.names = F,col.names = F)
write.table(x=Scer_regular_gene_expr_pca$x,file = paste0(output_workspace,"mtx_Scer_regular_gene_expr_pca.tsv"),sep = ",",na = "NA",row.names = F,col.names = F)

#Expression variance partitioned across the accessory and core genes 
v_comps <- c("Accessory Genome","Core Genome")
lst_comps <- list("comp_a"= c(v_comps[1],v_comps[1]),
                  "comp_b"=c(v_comps[1],v_comps[2]),
                  "comp_c"=c(v_comps[2],v_comps[2]),
                  "residual_comp"=c("Residual","Residual"))
#Across S cerevisiae strains (short timescale of evolution)
#import the variance partitioning bootstrap estimates 
df_Scer_Expression_varpart_BOOTSTRAP_estimates <- read.delim(paste0(output_workspace,"Table_Whole_Transcriptome_variance_explained_in_the_dataset_Bootstrap_Scerevisiae_1000p_strains.tsv"))
#weight variance explaines by PC importance
df_Scer_Expression_varpart_BOOTSTRAP_estimates <- df_Scer_Expression_varpart_BOOTSTRAP_estimates %>%
  pivot_longer(cols = (!id_subsample) & (!id_pc), names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(id_subsample,Component)  %>% 
  summarise(proportion_of_Expression_var_explained = sum(proportion_of_Expression_var_explained*df_Scer_regular_gene_expr_pc_eigenvalues$explained_variance[id_pc+1]))
#pivot to wide format before removing bootstrap without REML convergence
df_Scer_Expression_varpart_BOOTSTRAP_estimates <- df_Scer_Expression_varpart_BOOTSTRAP_estimates %>%
  pivot_wider(names_from = Component, values_from = proportion_of_Expression_var_explained)
#filter out the bootstrap without REML convergence
df_Scer_Expression_varpart_BOOTSTRAP_estimates <- subset(df_Scer_Expression_varpart_BOOTSTRAP_estimates,(comp_Accessory_genes_expression<=(1-residual_comp))&(comp_Core_genes_expression<=(1-residual_comp))&(residual_comp>=0))
#summarize bootstrap results
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- df_Scer_Expression_varpart_BOOTSTRAP_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(Component) %>% 
  summarise(mean_variance_explained = mean(proportion_of_Expression_var_explained), median_variance_explained = median(proportion_of_Expression_var_explained), se_variance_explained = sd(proportion_of_Expression_var_explained), IQR_variance_explained = IQR(proportion_of_Expression_var_explained))
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- vapply(X = df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- ifelse(test = df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="NA",yes = NA,no=df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start)
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- vapply(X = df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- ifelse(test = df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end=="NA",yes = NA,no=df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end)
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp <- vapply(X = 1:nrow(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) paste0(unique(c(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[the_i],df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp%in%c("NA","Residual")] <- NA
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination <- df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination[df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="Residual"] <- "Residual"
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$y_min <- vapply(X = 1:nrow(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) max(0,df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]-(1*df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$y_max <- vapply(X = 1:nrow(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) min(1,df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]+(1*df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2 <- NaN
df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2[df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- v_comps[2]
#dataframe for bottom part of the upset plot
df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats#subset(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))
for (i in 1:nrow(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats)){
  if ((!is.na(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]))&(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]!=df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[i])){
    v_current_comps <- strsplit(x = df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination[i],split = "&",fixed = T)[[1]]
    if ((length(v_current_comps)<=1)||(length(v_current_comps)>3)){
      next()
    }else if (length(v_current_comps)==2){
      #add end -> start
      current_subset1 <- df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset1$edge_end <- unname(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1)
    }else if (length(v_current_comps)==3){
      #add end -> start and middle -> end
      current_subset1 <- df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end2)
      current_subset2 <- df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset2$edge_start <- unname(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset2$edge_end <- unname(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(rbind(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1),current_subset2)
    }
  }else{
    next()
  }
}

#Upset plot Bootstrap
#upper part
upper_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_col(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained),fill="grey60") +
  geom_pointrange(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained,ymin = y_min,ymax = y_max ))+
  ylab("Proportion of gene expression variance explained") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=16),axis.text.y = element_text(size=11))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
#bottom part
bottom_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_bottom_upset_plot_Scer_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination)) + 
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_end,levels = c(v_comps,"Residual")),group = combination))+
  geom_line(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination))+
  xlab("Components") + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_text(size=16),axis.text.x = element_text(size=9))
#save upper and bottom parts of the upset plot
#ggsave(plot = upper_part_upsetPlot_Bootstrap_Varpart_gg, filename = "Scer_Expression_Varpart_Bootstrap_500_13_upsetPlot.svg", path=output_workspace, width = 20/2.54, height = 15/2.54,device = svg)
#ggsave(plot = bottom_part_upsetPlot_Bootstrap_Varpart_gg,filename = "Scer_bottom_part_upsetPlot_Bootstrap_500_13_gg.svg", path=output_workspace, width = 18/2.54, height = 8/2.54,device = svg)

#import the FULL dataset estimates
df_Scer_Expression_varpart_FULL_dataset_estimates <- read.delim(paste0(output_workspace,"Table_Whole_Transcriptome_variance_explained_in_the_dataset_ALL_Scerevisiae_1000p_strains.tsv"))
#weight variance explaines by PC importance
df_Scer_Expression_varpart_FULL_dataset_estimates_long <- df_Scer_Expression_varpart_FULL_dataset_estimates %>%
  pivot_longer(cols = (!id_subsample) & (!id_pc), names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(id_subsample,Component)  %>% 
  summarise(proportion_of_Expression_var_explained = sum(proportion_of_Expression_var_explained*df_Scer_regular_gene_expr_pc_eigenvalues$explained_variance[id_pc+1]))
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start <- vapply(X = df_Scer_Expression_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start <- ifelse(test = df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start=="NA",yes = NA,no=df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start)
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end <- vapply(X = df_Scer_Expression_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end <- ifelse(test = df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end=="NA",yes = NA,no=df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end)
df_Scer_Expression_varpart_FULL_dataset_estimates_long$name_comp <- vapply(X = 1:nrow(df_Scer_Expression_varpart_FULL_dataset_estimates_long),function(the_i) paste0(unique(c(df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start[the_i],df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Scer_Expression_varpart_FULL_dataset_estimates_long$name_comp[df_Scer_Expression_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Scer_Expression_varpart_FULL_dataset_estimates_long$name_comp[df_Scer_Expression_varpart_FULL_dataset_estimates_long$name_comp%in%c("NA","Residual")] <- NA
df_Scer_Expression_varpart_FULL_dataset_estimates_long$combination <- df_Scer_Expression_varpart_FULL_dataset_estimates_long$name_comp
df_Scer_Expression_varpart_FULL_dataset_estimates_long$combination[df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_start=="Residual"] <- "Residual"
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end2 <- NaN
df_Scer_Expression_varpart_FULL_dataset_estimates_long$edge_end2[df_Scer_Expression_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- v_comps[2]

#Euler diagram
eulerr_components <- subset(df_Scer_Expression_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$proportion_of_Expression_var_explained/df_Scer_regular_gene_expr_pc_eigenvalues$explained_variance_cum[nb_pcs_recapitulating_Scer_regular_gene_expr_pca]
names(eulerr_components) <- subset(df_Scer_Expression_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$combination
#Replace any potentially negative R^2 values with 0 for plotting consistency
eulerr_components[eulerr_components < 0] <- 0
#get %variance explained by dividing by sum(eulerr_components) as we selected the PCs exlaining AT LEAST 90% of transcription variance (and exactly sum(eulerr_components))
eulerr_components <- eulerr_components/sum(eulerr_components)

#Fit the Euler diagram to the data
fit <- eulerr::euler(eulerr_components)

svg(filename = paste0(output_workspace,"Scer_Expression_Varpart_FULL_dataset_1_969_Euler_Diagram.svg"),width = 8.8/2.54,height = 6/2.54)
#Plot the diagram
plot(fit,
     #Set the names of the three predictors (circles)
     labels = c("Accessory Genome","Core Genome", "Residual"),
     
     #Use the R^2 values as labels for the regions (rounded to 2 decimal places)
     quantities = lapply(eulerr_components, function(x) round(x, 2)),
     
     main = "Partitioning of S. cerevisiae's strains transcriptome variance (n=969)",
     
     #Customize colors and style
     fills = list(fill = c("#3b82f6", "#ef4444", "#22c55e"), alpha = 0.7),
     edges = list(lwd = 1.5, col = "black")
)
dev.off()

#Across Yeast species (long timescale of evolution)
#import the variance partitioning bootstrap estimates 
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates <- read.delim(paste0(output_workspace,"Table_Whole_Transcriptome_variance_explained_in_the_dataset_Bootstrap_Yeast_subphylum_subset.tsv"))
#weight variance explaines by PC importance
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates <- df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates %>%
  pivot_longer(cols = (!id_subsample) & (!id_pc), names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(id_subsample,Component)  %>% 
  summarise(proportion_of_Expression_var_explained = sum(proportion_of_Expression_var_explained*df_species_aggregate_gene_expr_pc_eigenvalues_YPD$explained_variance[id_pc+1]))
#pivot to wide format before removing bootstrap without REML convergence
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates <- df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates %>%
  pivot_wider(names_from = Component, values_from = proportion_of_Expression_var_explained)
#filter out the bootstrap without REML convergence
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates <- subset(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates,(comp_Accessory_genes_expression<=(1-residual_comp))&(comp_Core_genes_expression<=(1-residual_comp))&(residual_comp>=0))
#summarize bootstrap results
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates%>%
  pivot_longer(cols = !id_subsample, names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(Component) %>% 
  summarise(mean_variance_explained = mean(proportion_of_Expression_var_explained), median_variance_explained = median(proportion_of_Expression_var_explained), se_variance_explained = sd(proportion_of_Expression_var_explained), IQR_variance_explained = IQR(proportion_of_Expression_var_explained))
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- vapply(X = df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start <- ifelse(test = df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="NA",yes = NA,no=df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start)
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- vapply(X = df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end <- ifelse(test = df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end=="NA",yes = NA,no=df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end)
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) paste0(unique(c(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[the_i],df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp[df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp%in%c("NA","Residual")] <- NA
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination <- df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$name_comp
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination[df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start=="Residual"] <- "Residual"
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$y_min <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) max(0,df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]-(1*df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$y_max <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats),function(the_i) min(1,df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$mean_variance_explained[the_i]+(1*df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$se_variance_explained[the_i])),FUN.VALUE = 0.0)
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2 <- NaN
df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end2[df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$Component=="comp_e"] <- v_comps[2]
#dataframe for bottom part of the upset plot
df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats#subset(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))
for (i in 1:nrow(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats)){
  if ((!is.na(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]))&(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_start[i]!=df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$edge_end[i])){
    v_current_comps <- strsplit(x = df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats$combination[i],split = "&",fixed = T)[[1]]
    if ((length(v_current_comps)<=1)||(length(v_current_comps)>3)){
      next()
    }else if (length(v_current_comps)==2){
      #add end -> start
      current_subset1 <- df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset1$edge_end <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1)
    }else if (length(v_current_comps)==3){
      #add end -> start and middle -> end
      current_subset1 <- df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset1$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end2)
      current_subset2 <- df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]
      current_subset2$edge_start <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_end)
      current_subset2$edge_end <- unname(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats[i,]$edge_start)
      df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats <- rbind(rbind(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats, current_subset1),current_subset2)
    }
  }else{
    next()
  }
}

#Upset plot Bootstrap
#upper part
upper_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_col(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained),fill="grey60") +
  geom_pointrange(mapping = aes(x=reorder(combination,-mean_variance_explained),y=mean_variance_explained,ymin = y_min,ymax = y_max ))+
  ylab("Proportion of gene expression variance explained") + 
  theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.title.y = element_text(size=16),axis.text.y = element_text(size=11))+
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
#bottom part
bottom_part_upsetPlot_Bootstrap_Varpart_gg <- ggplot(data = subset(df_bottom_upset_plot_Yeast_subphylum_subset_Expression_varpart_BOOTSTRAP_estimates_summary_stats,!is.na(edge_start))) +
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination)) + 
  geom_point(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_end,levels = c(v_comps,"Residual")),group = combination))+
  geom_line(mapping = aes(x=reorder(gsub(pattern = "&",replacement="\n",x =combination),-mean_variance_explained),y=factor(edge_start,levels = c(v_comps,"Residual")),group = combination))+
  xlab("Components") + 
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x = element_text(size=16),axis.text.x = element_text(size=9))
#save upper and bottom parts of the upset plot
#ggsave(plot = upper_part_upsetPlot_Bootstrap_Varpart_gg, filename = "Yeast_subphylum_subset_Expression_Varpart_Bootstrap_500_13_upsetPlot.svg", path=output_workspace, width = 20/2.54, height = 15/2.54,device = svg)
#ggsave(plot = bottom_part_upsetPlot_Bootstrap_Varpart_gg,filename = "Yeast_subphylum_subset_bottom_part_upsetPlot_Bootstrap_500_13_gg.svg", path=output_workspace, width = 18/2.54, height = 8/2.54,device = svg)

#import the FULL dataset estimates
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates <- read.delim(paste0(output_workspace,"Table_Whole_Transcriptome_variance_explained_in_the_dataset_ALL_Yeast_subphylum_subset.tsv"))
#weight variance explaines by PC importance
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long <- df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates %>%
  pivot_longer(cols = (!id_subsample) & (!id_pc), names_to = "Component",values_to = "proportion_of_Expression_var_explained") %>%
  group_by(id_subsample,Component)  %>% 
  summarise(proportion_of_Expression_var_explained = sum(proportion_of_Expression_var_explained*df_species_aggregate_gene_expr_pc_eigenvalues_YPD$explained_variance[id_pc+1]))
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start <- vapply(X = df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][1],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start <- ifelse(test = df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start=="NA",yes = NA,no=df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start)
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end <- vapply(X = df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$Component,FUN = function(the_comp) ifelse(test = the_comp%in%names(lst_comps),yes = lst_comps[[the_comp]][2],no = "NA"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end <- ifelse(test = df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end=="NA",yes = NA,no=df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end)
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$name_comp <- vapply(X = 1:nrow(df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long),function(the_i) paste0(unique(c(df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start[the_i],df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end[the_i])),collapse = "&"),FUN.VALUE = "")
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$name_comp[df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- paste0(v_comps,collapse = "&")
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$name_comp[df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$name_comp%in%c("NA","Residual")] <- NA
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$combination <- df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$name_comp
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$combination[df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_start=="Residual"] <- "Residual"
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end2 <- NaN
df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$edge_end2[df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long$Component=="comp_e"] <- v_comps[2]

#Euler diagram
eulerr_components <- subset(df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$proportion_of_Expression_var_explained
names(eulerr_components) <- subset(df_Yeast_subphylum_subset_Expression_varpart_FULL_dataset_estimates_long,!is.na(edge_start))$combination
#Replace any potentially negative R^2 values with 0 for plotting consistency
eulerr_components[eulerr_components < 0] <- 0
#get %variance explained by dividing by sum(eulerr_components) as we selected the PCs exlaining AT LEAST 90% of transcription variance (and exactly sum(eulerr_components))
eulerr_components <- eulerr_components/sum(eulerr_components)

#Fit the Euler diagram to the data
fit <- eulerr::euler(eulerr_components)

svg(filename = paste0(output_workspace,"Yeast_subphylum_subset_Expression_Varpart_FULL_dataset_1_16_Euler_Diagram.svg"),width = 8.8/2.54,height = 6/2.54)
#Plot the diagram
plot(fit,
     #Set the names of the three predictors (circles)
     labels = c("Accessory Genome","Core Genome", "Residual"),
     
     #Use the R^2 values as labels for the regions (rounded to 2 decimal places)
     quantities = lapply(eulerr_components, function(x) round(x, 2)),
     
     main = "Partitioning of Yeast species transcriptome variance (n=16)",
     
     #Customize colors and style
     fills = list(fill = c("#3b82f6", "#ef4444", "#22c55e"), alpha = 0.7),
     edges = list(lwd = 1.5, col = "black")
)
dev.off()


#**********************************************************************Whole-Transcriptome evolution analysis*********************************************************************#
#Whole-Transcriptome evolution analysis
  #extract distances (time = phylo distance)
v_Transcriptome_similarity_of_almost_21_spcs_YPD <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairwise_expression_similarity_YPD)
v_phylo_dist_of_almost_21_spcs_YPD <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_phylo_dist_the_21_sps_YPD)
v_most_influencial_Orthogroup_expression_similarity <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(1/vegdist(x = as.matrix(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,6543]),method = "euclidean")))
v_most_influencial_Orthogroup_expression_similarity[!is.finite(v_most_influencial_Orthogroup_expression_similarity)] <- NA
v_fitness_similarity_of_almost_21_spcs_YPD <- 1/v_dist_fitness_21_spcs_YPD
v_fitness_similarity_of_almost_21_spcs_YPD[!is.finite(v_fitness_similarity_of_almost_21_spcs_YPD)] <- NA
mtx_pairs_of_almost_21_spcs_YPD <- matrix(NA,nrow = nrow(mtx_phylo_dist_the_21_sps_YPD),ncol = ncol(mtx_phylo_dist_the_21_sps_YPD))
rownames(mtx_pairs_of_almost_21_spcs_YPD) <- gsub(pattern = "_",replacement = " ",rownames(mtx_phylo_dist_the_21_sps_YPD))
colnames(mtx_pairs_of_almost_21_spcs_YPD) <- gsub(pattern = "_",replacement = " ",colnames(mtx_phylo_dist_the_21_sps_YPD))
for (i in 1:nrow(mtx_phylo_dist_the_21_sps_YPD)){
  for (j in 1:ncol(mtx_phylo_dist_the_21_sps_YPD)){
    mtx_pairs_of_almost_21_spcs_YPD[i,j] <- paste0(gsub(pattern = "_",replacement = " ",rownames(mtx_phylo_dist_the_21_sps_YPD)[i]),";",gsub(pattern = "_",replacement = " ",colnames(mtx_phylo_dist_the_21_sps_YPD)[j]))
  }
}
df_YPD_almost_21_spcs_pairs_features_comparison <- data.frame(sp1=unname(vapply(X = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairs_of_almost_21_spcs_YPD),FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,";"))[1],FUN.VALUE = "")),
                                                              sp2=unname(vapply(X = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairs_of_almost_21_spcs_YPD),FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,";"))[2],FUN.VALUE = "")),
                                                              phylo_dist=v_phylo_dist_of_almost_21_spcs_YPD,
                                                              Transcriptome_similarity=v_Transcriptome_similarity_of_almost_21_spcs_YPD,
                                                              Transcriptome_similarity_based_on_present_genes_only = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairwise_expression_similarity_based_on_present_genes_YPD),
                                                              Most_influencial_Orthogroup_expression_similarity=v_most_influencial_Orthogroup_expression_similarity,
                                                              Fitness_similarity=v_fitness_similarity_of_almost_21_spcs_YPD) 
df_YPD_almost_21_spcs_pairs_features_comparison$label_species <- paste0(df_YPD_almost_21_spcs_pairs_features_comparison$sp1," vs ",df_YPD_almost_21_spcs_pairs_features_comparison$sp2)
df_YPD_almost_21_spcs_pairs_features_comparison <- df_YPD_almost_21_spcs_pairs_features_comparison %>%
  dplyr::arrange((phylo_dist))
  #define species clade / taxonomic order
v_order_almost_21_spcs_YPD <- c("Pichiales", "Dipodascales", "Serinales", "Serinales",
                                "Phaffomycetales", "Phaffomycetales", "Phaffomycetales", "Saccharomycodales",
                                "Saccharomycodales", "Serinales", "Serinales", "Alaninales",
                                "Alaninales", "Alloascoideales", "Phaffomycetales", "Phaffomycetales")
names(v_order_almost_21_spcs_YPD) <- rownames(mtx_pairs_of_almost_21_spcs_YPD)
df_YPD_almost_21_spcs_pairs_features_comparison$clade_order_sp1 <- v_order_almost_21_spcs_YPD[df_YPD_almost_21_spcs_pairs_features_comparison$sp1]
df_YPD_almost_21_spcs_pairs_features_comparison$clade_order_sp2 <- v_order_almost_21_spcs_YPD[df_YPD_almost_21_spcs_pairs_features_comparison$sp2]
    #clade color
v_clade_color <- c("#6A3906","#008000","#FF2800","#00E4FF","#FF8200","#12958A","#FFD200")
names(v_clade_color) <- unique(sort(v_order_almost_21_spcs_YPD))
    #create a color palette to identify species
my_palette <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A",
                "#FF7F00", "black", "gold1", "skyblue2",
                "palegreen2", "#FDBF6F", "gray70", "maroon",
                "orchid1", "darkturquoise", "darkorange4", "brown")
names(my_palette) <- rownames(mtx_pairs_of_almost_21_spcs_YPD)

  #function to get the equation and R^2 of a linear model as a string
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == b %.% italic(x) + a*","~~italic(R)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

  #Transcriptome similarity vs Phylogenetic distance
    #color plot based on species pairs (fill = sp1, color = sp2)
Transcriptome_similarity_vs_divergence_time_colored_by_sps_pairs_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison, aes(x = phylo_dist, y = Transcriptome_similarity, color = sp1, fill = sp2)) +
  geom_point(shape=21,size=3,stroke=2) + 
  scale_fill_manual(values=my_palette) + 
  scale_colour_manual(values=my_palette)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") +
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Species") +
  guides(color = "none")

    #color plot based on species orders (fill = clade_order_sp1, color = clade_order_sp2)
Transcriptome_similarity_vs_divergence_time_colored_by_sps_orders_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison, aes(x = phylo_dist, y = Transcriptome_similarity, col = clade_order_sp1, fill = clade_order_sp2)) +
  geom_point(shape=21,size=3,stroke=2) + 
  scale_fill_manual(values=v_clade_color) + 
  scale_colour_manual(values=v_clade_color)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Order") +
  guides(color = "none") + 
  theme(legend.key.size = unit(1, "cm")) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
Transcriptome_similarity_vs_divergence_time_gg <- grid.arrange(Transcriptome_similarity_vs_divergence_time_colored_by_sps_pairs_gg,Transcriptome_similarity_vs_divergence_time_colored_by_sps_orders_gg,ncol=1)
ggsave(plot = Transcriptome_similarity_vs_divergence_time_gg,filename = "Yeast_subphylum_Transcriptome_similarity_vs_divergence_time_gg.svg", path=output_workspace, width = 25/2.54, height = 25/2.54,device = svg)
ggsave(plot = Transcriptome_similarity_vs_divergence_time_gg,filename = "Yeast_subphylum_Transcriptome_similarity_vs_divergence_time_gg.png", path=output_workspace, width = 25/2.54, height = 25/2.54,device = png,dpi = 600)
    #linear model summary
summary(lm(data = df_YPD_almost_21_spcs_pairs_features_comparison, formula = Transcriptome_similarity~phylo_dist))

    #re-plot but with regression line
df_YPD_almost_21_spcs_pairs_features_comparison$x <- df_YPD_almost_21_spcs_pairs_features_comparison$phylo_dist
df_YPD_almost_21_spcs_pairs_features_comparison$y <- df_YPD_almost_21_spcs_pairs_features_comparison$Transcriptome_similarity
lbl_eqn_rsq_transc_sim_vs_phylo_dist <- lm_eqn(df = df_YPD_almost_21_spcs_pairs_features_comparison)
Transcriptome_similarity_vs_divergence_time_colored_by_sps_orders_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison) +
  geom_point(mapping = aes(x = phylo_dist, y = Transcriptome_similarity, col = clade_order_sp1, fill = clade_order_sp2),shape=21,size=3,stroke=2) + 
  geom_smooth(mapping = aes(x = phylo_dist, y = Transcriptome_similarity),method = "lm", se=FALSE, color="black",lty=2) +
  geom_text(x = 3, y = 0.8, label = lbl_eqn_rsq_transc_sim_vs_phylo_dist, parse = TRUE, cex = 6) + 
  scale_fill_manual(values=v_clade_color) + 
  scale_colour_manual(values=v_clade_color)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Order",title="All orthogroups included") + #absent orthogroups have an expression of 0 to illustrate the impact of gene content variation and regulatory SNPs on expression
  guides(color = "none") + 
  theme(legend.key.size = unit(1, "cm")) +
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))

#Transcriptome similarity (ONLY BASED ON GENES DETECTED IN THE GENOME) vs Phylogenetic distance
  #color plot based on species pairs (fill = sp1, color = sp2)
Transcriptome_similarity_based_on_present_genes_only_vs_divergence_time_colored_by_sps_pairs_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison, aes(x = phylo_dist, y = Transcriptome_similarity_based_on_present_genes_only)) +
  geom_point(shape=21,size=3,stroke=2) + 
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") +
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Species") +
  guides(color = "none")

#color plot based on species orders (fill = clade_order_sp1, color = clade_order_sp2)
df_YPD_almost_21_spcs_pairs_features_comparison$x <- df_YPD_almost_21_spcs_pairs_features_comparison$phylo_dist
df_YPD_almost_21_spcs_pairs_features_comparison$y <- df_YPD_almost_21_spcs_pairs_features_comparison$Transcriptome_similarity_based_on_present_genes_only
lbl_eqn_rsq_transc_sim_vs_phylo_dist <- lm_eqn(df = df_YPD_almost_21_spcs_pairs_features_comparison)
Transcriptome_similarity_based_on_present_genes_only_vs_divergence_time_colored_by_sps_orders_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison) +
  geom_point(mapping = aes(x = phylo_dist, y = Transcriptome_similarity_based_on_present_genes_only),shape=21,size=3,stroke=2) + 
  geom_smooth(mapping = aes(x = phylo_dist, y = Transcriptome_similarity_based_on_present_genes_only),method = "lm", se=FALSE, color="black",lty=2) +
  geom_text(x = 3, y = 0.8, label = lbl_eqn_rsq_transc_sim_vs_phylo_dist, parse = TRUE, cex = 6) + 
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Order",title = "Only orthogroups that are in both genomes are included") + #absent orthgroups are not accounted for in the similarity computation to highlight the impact of regulatory SNPs only
  guides(color = "none") + 
  theme(legend.key.size = unit(1, "cm")) + 
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
df_YPD_almost_21_spcs_pairs_features_comparison$x <- NULL
df_YPD_almost_21_spcs_pairs_features_comparison$y <- NULL

Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg <- grid.arrange(Transcriptome_similarity_vs_divergence_time_colored_by_sps_orders_gg,Transcriptome_similarity_based_on_present_genes_only_vs_divergence_time_colored_by_sps_orders_gg,ncol=1)
ggsave(plot = Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg,filename = "Yeast_subphylum_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg.svg", path=output_workspace, width = 25/2.54, height = 25/2.54,device = svg)
ggsave(plot = Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg,filename = "Yeast_subphylum_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg.png", path=output_workspace, width = 25/2.54, height = 25/2.54,device = png,dpi = 600)
#linear model summary
summary(lm(data = df_YPD_almost_21_spcs_pairs_features_comparison, formula = Transcriptome_similarity_based_on_present_genes_only~phylo_dist))

#    #interactive plot
#plot_ly(data = df_YPD_almost_21_spcs_pairs_features_comparison, x = ~phylo_dist, y = ~Transcriptome_similarity, text = ~label_species,
#        type = 'scatter', mode = 'markers', marker = list(
#          line = list(color = ~v_clade_color[clade_order_sp1], width = 3), #Point outline color (CONTOUR)
#          color = ~v_clade_color[clade_order_sp2], #Point FILL color
#          size = 10
#        )) %>%
#  layout(title = "Yeast subphylum transcriptome similarity in function of divergence time and colored by clades",
#         xaxis = list(title = "Phylogenetic similarity (divergence time or #substitutions per site)"),
#         yaxis = list(title = "Transcriptome similarity (R^2)"))

  #Most influential orthogroup (OG0007735) expression vs phylo distance (divergence time)
    #color plot based on clade pairs (fill = sp1, color = sp2)
Similarity_in_the_most_influential_transcript_vs_divergence_time_colored_by_sps_pairs_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison, aes(x = phylo_dist, y = Most_influencial_Orthogroup_expression_similarity, col = sp1, fill = sp2)) +
  geom_point(shape=21,size=3,stroke=2) + 
  scale_fill_manual(values=my_palette) + 
  scale_colour_manual(values=my_palette)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Similarity in the transcript level of the most influential gene\n(OG0007735:Cell wall oganization protein Zeo1)'~(R^2)))+
  labs(fill = "Species") +
  guides(color = "none")
    #color plot based on species orders (fill = clade_order_sp1, color = clade_order_sp2)
Similarity_in_the_most_influential_transcript_vs_divergence_time_colored_by_sps_orders_gg <- ggplot(data = df_YPD_almost_21_spcs_pairs_features_comparison, aes(x = phylo_dist, y = Most_influencial_Orthogroup_expression_similarity, col = clade_order_sp1, fill = clade_order_sp2)) +
  geom_point(shape=21,size=3,stroke=2) + 
  scale_fill_manual(values=v_clade_color) + 
  scale_colour_manual(values=v_clade_color)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Similarity in the transcript level of the most influential gene\n(OG0007735:Cell wall oganization protein Zeo1)'~(R^2)))+
  labs(fill = "Species") +
  guides(color = "none")
Similarity_in_the_most_influential_transcript_vs_divergence_time_gg <- grid.arrange(Similarity_in_the_most_influential_transcript_vs_divergence_time_colored_by_sps_pairs_gg,Similarity_in_the_most_influential_transcript_vs_divergence_time_colored_by_sps_orders_gg,ncol=1)
ggsave(plot = Similarity_in_the_most_influential_transcript_vs_divergence_time_gg,filename = "Yeast_subphylum_Similarity_in_the_most_influential_transcript_vs_divergence_time_gg.svg", path=output_workspace, width = 25/2.54, height = 25/2.54,device = svg)
#    #interactive plot
#plot_ly(data = df_YPD_almost_21_spcs_pairs_features_comparison, x = ~phylo_dist, y = ~Most_influencial_Orthogroup_expression_similarity, text = ~label_species,
#        type = 'scatter', mode = 'markers', marker = list(
#          line = list(color = ~v_clade_color[clade_order_sp1], width = 3), #Point outline color (CONTOUR)
#          color = ~v_clade_color[clade_order_sp2], #Point FILL color
#          size = 10
#        )) %>%
#  layout(title = "Similarity in the transcript level of the most influential gene in function of divergence time and colored by clades",
#         xaxis = list(title = "Phylogenetic similarity (divergence time or #substitutions per site)"),
#         yaxis = list(title = "Similarity in the transcript level of the most influential gene (R^2; OG0007735:Cell wall oganization protein Zeo1)"))

  #log10(Fitness similarity) vs Transcriptome similarity + Phylogenetic distance is a small correlation (random effect need to be controlled for?)
summary(lm(data = df_YPD_almost_21_spcs_pairs_features_comparison, formula = log10(Fitness_similarity)~Transcriptome_similarity+phylo_dist))

  #For Scerevisiae strains
v_Transcriptome_similarity_of_Scer_strains <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(cor(t(mtx_Scer_expression_all_genes))^2)[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect])
v_phylo_dist_of_Scer_strains <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_phylo_dist_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect])
v_fitness_similarity_of_Scer_strains <- 1/v_dist_Scer_strains_fitness
v_fitness_similarity_of_Scer_strains[!is.finite(v_fitness_similarity_of_Scer_strains)] <- NA
mtx_pairs_of_Scer_strains <- matrix(NA,nrow = length(v_Scer_matrices_strain_intersect),ncol = length(v_Scer_matrices_strain_intersect))
for (i in 1:length(v_Scer_matrices_strain_intersect)){
  for (j in 1:length(v_Scer_matrices_strain_intersect)){
    mtx_pairs_of_Scer_strains[i,j] <- paste0(v_Scer_matrices_strain_intersect[i],";",v_Scer_matrices_strain_intersect[j])
  }
}
df_Scer_strains_pairs_features_comparison <- data.frame(Scer_strain1=unname(vapply(X = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairs_of_Scer_strains),FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,";"))[1],FUN.VALUE = "")),
                                                              Scer_strain2=unname(vapply(X = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairs_of_Scer_strains),FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,";"))[2],FUN.VALUE = "")),
                                                              phylo_dist=v_phylo_dist_of_Scer_strains,
                                                              Transcriptome_similarity=v_Transcriptome_similarity_of_Scer_strains,
                                                              Transcriptome_similarity_based_on_present_genes_only = get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairwise_expression_similarity_based_on_present_genes_in_Scer_strains[v_Scer_matrices_strain_intersect,v_Scer_matrices_strain_intersect]),
                                                              Fitness_similarity=v_fitness_similarity_of_Scer_strains) 
df_Scer_strains_pairs_features_comparison$label_species <- paste0(df_Scer_strains_pairs_features_comparison$Scer_strain1," vs ",df_Scer_strains_pairs_features_comparison$Scer_strain2)
df_Scer_strains_pairs_features_comparison <- df_Scer_strains_pairs_features_comparison %>%
  dplyr::arrange((phylo_dist))


#Scer strains transcriptome similarity based on all genes
df_Scer_strains_pairs_features_comparison$x <- df_Scer_strains_pairs_features_comparison$phylo_dist
df_Scer_strains_pairs_features_comparison$y <- df_Scer_strains_pairs_features_comparison$Transcriptome_similarity
lbl_eqn_rsq_transc_sim_vs_phylo_dist <- lm_eqn(df = df_Scer_strains_pairs_features_comparison)
Scer_Transcriptome_similarity_vs_divergence_time_gg <- ggplot(data = df_Scer_strains_pairs_features_comparison) +
  geom_point(mapping = aes(x = phylo_dist, y = Transcriptome_similarity, col = Scer_strain1, fill = Scer_strain2),shape=21,size=3,stroke=2) + 
  geom_smooth(mapping = aes(x = phylo_dist, y = Transcriptome_similarity),method = "lm", se=FALSE, color="black",lty=2) +
  geom_text(x = 3, y = 0.8, label = lbl_eqn_rsq_transc_sim_vs_phylo_dist, parse = TRUE, cex = 6) + 
  scale_fill_manual(values=v_clade_color) + 
  scale_colour_manual(values=v_clade_color)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Strain2",title = "Only orthogroups that are in both genomes are included") + #absent orthgroups are not accounted for in the similarity computation to highlight the impact of regulatory SNPs only
  guides(color = "none") + 
  theme(legend.key.size = unit(1, "cm")) + 
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
df_Scer_strains_pairs_features_comparison$x <- NULL
df_Scer_strains_pairs_features_comparison$y <- NULL

  #linear model summary
summary(lm(data = df_Scer_strains_pairs_features_comparison, formula = Transcriptome_similarity~phylo_dist))

#Scer strains transcriptome similarity based on genes shared by the strain pair
df_Scer_strains_pairs_features_comparison$x <- df_Scer_strains_pairs_features_comparison$phylo_dist
df_Scer_strains_pairs_features_comparison$y <- df_Scer_strains_pairs_features_comparison$Transcriptome_similarity_based_on_present_genes_only
lbl_eqn_rsq_transc_sim_vs_phylo_dist <- lm_eqn(df = df_Scer_strains_pairs_features_comparison)
Scer_Transcriptome_similarity_based_on_present_genes_only_vs_divergence_time_gg <- ggplot(data = df_Scer_strains_pairs_features_comparison) +
  geom_point(mapping = aes(x = phylo_dist, y = Transcriptome_similarity_based_on_present_genes_only, col = Scer_strain1, fill = Scer_strain2),shape=21,size=3,stroke=2) + 
  geom_smooth(mapping = aes(x = phylo_dist, y = Transcriptome_similarity_based_on_present_genes_only),method = "lm", se=FALSE, color="black",lty=2) +
  geom_text(x = 3, y = 0.8, label = lbl_eqn_rsq_transc_sim_vs_phylo_dist, parse = TRUE, cex = 6) + 
  scale_fill_manual(values=v_clade_color) + 
  scale_colour_manual(values=v_clade_color)+
  xlab("Phylogenetic distance (divergence time or #substitutions per site)") + 
  ylab(bquote('Transcriptome similarity'~(R^2)))+
  labs(fill = "Strain2",title = "Only orthogroups that are in both genomes are included") + #absent orthgroups are not accounted for in the similarity computation to highlight the impact of regulatory SNPs only
  guides(color = "none") + 
  theme(legend.key.size = unit(1, "cm")) + 
  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
df_Scer_strains_pairs_features_comparison$x <- NULL
df_Scer_strains_pairs_features_comparison$y <- NULL
  #linear model summary
summary(lm(data = df_Scer_strains_pairs_features_comparison, formula = Transcriptome_similarity_based_on_present_genes_only~phylo_dist))

#Grid figure with both panels
Scer_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg <- grid.arrange(Scer_Transcriptome_similarity_vs_divergence_time_gg,Scer_Transcriptome_similarity_based_on_present_genes_only_vs_divergence_time_gg,ncol=1)
ggsave(plot = Scer_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg,filename = "Yeast_subphylum_Scer_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg.svg", path=output_workspace, width = 25/2.54, height = 25/2.54,device = svg)
ggsave(plot = Scer_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg,filename = "Yeast_subphylum_Scer_Transcriptome_similarity_based_either_on_all_genes_or_only_based_on_present_genes_vs_divergence_time_gg.png", path=output_workspace, width = 25/2.54, height = 25/2.54,device = png,dpi = 600)

#******************************************************************************Per-Kegg pathway analysis + look for KEGG class enrichment********************************************************#
#Per-Kegg pathway analysis + look for KEGG class enrichment (CONSIDER THAT a single orthogroup / ko_id can be associated to multiple pathways and that each pathway has multiple orthogroups / ko_ids). THEREFORE, USE df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY
  #Kegg pathway orthogroups correlation between all pair of species recovered from YPD
    #correspondence between ALL possible Ko pathways and pathway titles
v_KO_pathway_to_pathway_title <- unique(df_KO_ID_to_KEGG_PATHWAY[,c("ko_pathway","pathway_title")])$pathway_title
names(v_KO_pathway_to_pathway_title) <- unique(df_KO_ID_to_KEGG_PATHWAY[,c("ko_pathway","pathway_title")])$ko_pathway
    #tip distances
v_tip_dists <- tip_distances[v_almost_21_spcs_to_their_names_with_underscore[colnames(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs)]]
    #Find ONLY the YPD KO pathways
v_lst_KO_pathways_YPD <- sort(unique(subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,condition=="YPD")$ko_pathway))

#find clusters of pathway size (#of shared orthogroups in species pairs) AND pathway size bins
mtx_KO_pathway_size_in_sp_pairs <- matrix(NA,nrow = length(v_lst_KO_pathways_YPD),ncol = length(v_lst_all_pairs_used_for_phylo_dist))
rownames(mtx_KO_pathway_size_in_sp_pairs) <- v_lst_KO_pathways_YPD
colnames(mtx_KO_pathway_size_in_sp_pairs) <- v_lst_all_pairs_used_for_phylo_dist
i <- 1
for (current_pathway in rownames(mtx_KO_pathway_size_in_sp_pairs)){
  mtx_KO_pathway_size_in_sp_pairs[current_pathway,] <- unname(vapply(X = v_lst_all_pairs_used_for_phylo_dist,FUN = function(current_sp_pair) length(intersect(subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,(!is.na(bool_is_detected_in_genome))&(bool_is_detected_in_genome)&(species_smpl_lbl==unlist(strsplit(current_sp_pair,";"))[1])&(ko_pathway==current_pathway))$Orthogroup, subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,(!is.na(bool_is_detected_in_genome))&(bool_is_detected_in_genome)&(species_smpl_lbl==unlist(strsplit(current_sp_pair,";"))[2])&(ko_pathway==current_pathway))$Orthogroup)),FUN.VALUE = 0))
  print(paste0(i," ko_pathway(s) done out of ",nrow(mtx_KO_pathway_size_in_sp_pairs),"!"))
  i <- i + 1
}


#draw 2-way clustered heatmap for PATHWAY SIZE without the dendrogramm, without the row/column labels, but with a side color bar showing the pathway class
  #filter out pathways that are never shared at the genomic level by species pairs for which we have expression data in YPD
v_lst_KO_pathways_shared_by_YPD_species_pairs <- rownames(mtx_KO_pathway_size_in_sp_pairs) [rowSums(mtx_KO_pathway_size_in_sp_pairs,na.rm = T)!=0]
mtx_KO_pathway_size_in_sp_pairs <- mtx_KO_pathway_size_in_sp_pairs[rowSums(mtx_KO_pathway_size_in_sp_pairs,na.rm = T)!=0,]
  #Bray-Curtis distance for rows
row_dist <- vegdist(mtx_KO_pathway_size_in_sp_pairs, method = "bray")
  #Bray-Curtis distance for columns (transpose the matrix first)
col_dist <- vegdist(t(mtx_KO_pathway_size_in_sp_pairs), method = "bray")
  #Ward clustering for rows
row_clust <- hclust(row_dist, method = "ward.D2") #Use ward.D2 for Ward's method
  #Ward clustering for columns
col_clust <- hclust(col_dist, method = "ward.D2")
  #get the HEATMAP clustered rownames
hmp_KO_pathway_size_in_sp_pairs <- heatmap.2(log10(mtx_KO_pathway_size_in_sp_pairs+1),
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust),
          Colv = as.dendrogram(col_clust),
          dendrogram = "none",
          trace = "none",
          main = NULL,
          xlab = "Species pair",
          ylab="KEGG Orthology Pathway",
          key.par = list(cex = 0.55),
          key.title = NA,
          key.xlab = "log10(pathway size+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          labRow=NA,
          labCol=NA,
          margins = c(1.5, 1.5),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2))
clustered_rownames <- rownames(mtx_KO_pathway_size_in_sp_pairs)[hmp_KO_pathway_size_in_sp_pairs$rowInd]
v_lst_ko_pathway_classes <- sort(unique(df_KO_ID_to_KEGG_PATHWAY$pathway_class))
palette_pathway_classes <- brewer.pal(n = length(v_lst_ko_pathway_classes),name = "Set1")
names(palette_pathway_classes) <- v_lst_ko_pathway_classes

mtx_pathway_to_pthw_class_association <- as.matrix(table(df_KO_ID_to_KEGG_PATHWAY$ko_pathway,df_KO_ID_to_KEGG_PATHWAY$pathway_class))
v_KO_pathway_to_most_frequent_pathway_class <- vapply(X = v_lst_KO_pathways_shared_by_YPD_species_pairs,FUN = function(current_pathway) colnames(mtx_pathway_to_pthw_class_association)[mtx_pathway_to_pthw_class_association[current_pathway,]==max(mtx_pathway_to_pthw_class_association[current_pathway,])],FUN.VALUE = "")

mtx_pathway_to_pthw_group_association <- as.matrix(table(df_KO_ID_to_KEGG_PATHWAY$ko_pathway,df_KO_ID_to_KEGG_PATHWAY$pathway_group))
v_KO_pathway_to_most_frequent_pathway_group <- vapply(X = v_lst_KO_pathways_shared_by_YPD_species_pairs,FUN = function(current_pathway) colnames(mtx_pathway_to_pthw_group_association)[mtx_pathway_to_pthw_group_association[current_pathway,]==max(mtx_pathway_to_pthw_group_association[current_pathway,])],FUN.VALUE = "")

  #draw the heatmap again but this time with the pathway classes in a side color bar
png(filename = paste0(output_workspace,"Heatmap_pathway_size_in_sp_pairs.png"),width = 25,height = 25,units = "cm",res = 300)
heatmap.2(log10(mtx_KO_pathway_size_in_sp_pairs+1),
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust),
          Colv = as.dendrogram(col_clust),
          dendrogram = "none",
          trace = "none",
          main = NULL,
          xlab = "Species pair",
          ylab="KEGG Orthology Pathway",
          key.par = list(cex = 0.55),
          key.title = NA,
          key.xlab = "log10(pathway size+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          labRow=NA,
          labCol=NA,
          margins = c(1.5, 1.5),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2)) #,RowSideColors = unname(palette_pathway_classes[v_KO_pathway_to_most_frequent_pathway_class[clustered_rownames]])
dev.off()
  #Draw silhouette score profile
dist_matrix <- row_dist #Calculate distance matrix
hclust_result <- row_clust #Perform hierarchical clustering
k_range <- 2:12
cluster_assignments <- list()
for (k in k_range) {
  cluster_assignments[[as.character(k)]] <- cutree(hclust_result, k = k)
}
silhouette_scores <- numeric(length(k_range))
for (i in seq_along(k_range)) {
  k <- k_range[i]
  sil_obj <- cluster::silhouette(cluster_assignments[[as.character(k)]], dist_matrix)
  silhouette_scores[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"Silhouette_scores_pathway_size_in_sp_pairs.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range, silhouette_scores, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
  #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_ko_pathways <- (2:12)[silhouette_scores==max(silhouette_scores)]
ko_pathways_clusters_vector <- cutree(row_clust, k = nb_clusters_ko_pathways)

  #View the row cluster assignments
print(ko_pathways_clusters_vector)

  #create pathway-to-bins correspondence vectors
v_pathway_to_size_clusters_or_bins <- ko_pathways_clusters_vector
names(v_pathway_to_size_clusters_or_bins) <- rownames(mtx_KO_pathway_size_in_sp_pairs)
  #Chisq test between Pathway classes and pathway size cluster
    #contingency table
contingency_table_chisq_association_pthw_class_and_size_cluster  <- table(unname(v_KO_pathway_to_most_frequent_pathway_class),unname(v_pathway_to_size_clusters_or_bins[names(v_KO_pathway_to_most_frequent_pathway_class)]))
contingency_table_chisq_association_pthw_class_and_size_cluster 
    #pvalue
pval_chisq_association_pthw_class_and_size_cluster <- chisq.test(contingency_table_chisq_association_pthw_class_and_size_cluster)$p.value
pval_chisq_association_pthw_class_and_size_cluster
    #Cramer's V^2
squared_cramerV_association_pthw_class_and_size_cluster <- lsr::cramersV(contingency_table_chisq_association_pthw_class_and_size_cluster)^2
squared_cramerV_association_pthw_class_and_size_cluster

#create a matrix of species pairs distances based on pathway-specific expression data
mtx_KO_pathway_specific_expression_distances <- matrix(NA,nrow = length(v_lst_KO_pathways_shared_by_YPD_species_pairs),ncol = length(v_lst_all_pairs_used_for_phylo_dist))
rownames(mtx_KO_pathway_specific_expression_distances) <- v_lst_KO_pathways_shared_by_YPD_species_pairs
colnames(mtx_KO_pathway_specific_expression_distances) <- v_lst_all_pairs_used_for_phylo_dist
for (current_pathway in rownames(mtx_KO_pathway_specific_expression_distances)){
  v_current_ko_pathway_Orthogroups <- unique(subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,ko_pathway==current_pathway)$Orthogroup)
  mtx_KO_pathway_specific_expression_distances[current_pathway,] <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_current_ko_pathway_Orthogroups], method = "bray")))
}
  #filter out KO pathways that are not shared by any species pairs
mtx_KO_pathway_specific_expression_distances <- mtx_KO_pathway_specific_expression_distances[rowSums(mtx_KO_pathway_specific_expression_distances,na.rm = T)!=0,]

#draw 2-way clustered heatmap for PATHWAY-SPECIFIC EXPRESSION DISTANCES IN SPECIES PAIR with a side color bar showing the pathway class
  #Bray-Curtis distance for rows
row_dist_sp_pair_expression_dist <- vegdist(mtx_KO_pathway_specific_expression_distances, method = "bray",na.rm=T)

  #Bray-Curtis distance for columns (transpose the matrix first)
col_dist_sp_pair_expression_dist <- vegdist(t(mtx_KO_pathway_specific_expression_distances), method = "bray",na.rm=T)

  #Ward clustering for rows
row_clust_sp_pair_expression_dist <- hclust(row_dist_sp_pair_expression_dist, method = "ward.D2") #Use ward.D2 for Ward's method

  #Ward clustering for columns
col_clust_sp_pair_expression_dist <- hclust(col_dist_sp_pair_expression_dist, method = "ward.D2")

  #get the HEATMAP clustered rownames
hmp_sp_pair_expression_dist <- heatmap.2(log10(mtx_KO_pathway_specific_expression_distances+1),
                                         distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
                                         hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
                                         Rowv = as.dendrogram(row_clust_sp_pair_expression_dist),
                                         Colv = as.dendrogram(col_clust_sp_pair_expression_dist),
                                         dendrogram = "none",
                                         trace = "none",
                                         main = NULL,
                                         xlab = "Species pair",
                                         ylab="KEGG Orthology Pathway",
                                         key.par = list(cex = 0.55),
                                         key.title = NA,
                                         key.xlab = "log10(Expression distance+1)",
                                         key.ylab = "Count",
                                         col= brewer.pal(n = 9, name = "YlOrRd"),
                                         labRow=NA,
                                         labCol=NA,
                                         margins = c(1.5, 1.5),
                                         lwid = c(1.2, 8),
                                         lhei = c(0.5, 2))
clustered_rownames_sp_pair_expression_dist <- rownames(mtx_KO_pathway_specific_expression_distances)[hmp_sp_pair_expression_dist$rowInd]
  #draw the heatmap again but this time with the pathway classes in a side color bar
png(filename = paste0(output_workspace,"Heatmap_sp_pair_expression_dist.png"),width = 25,height = 25,units = "cm",res = 300)
heatmap.2(log10(mtx_KO_pathway_specific_expression_distances+1),
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust_sp_pair_expression_dist),
          Colv = as.dendrogram(col_clust_sp_pair_expression_dist),
          dendrogram = "row",
          trace = "none",
          main = NULL,
          xlab = "Species pair",
          ylab="KEGG Orthology Pathway",
          key.par = list(cex = 0.55),
          key.title = NA,
          key.xlab = "log10(Expression distance+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          labRow=NA,
          labCol=NA,
          margins = c(1.5, 1.5),
          lwid = c(1.2, 5),
          lhei = c(0.5, 2),
          RowSideColors = unname(palette_pathway_classes[v_KO_pathway_to_most_frequent_pathway_class[clustered_rownames_sp_pair_expression_dist]])) #
dev.off()
  #Draw silhouette score profile
k_range_sp_pair_expression_dist <- seq(2,200,5)
cluster_assignments_sp_pair_expression_dist <- list()
for (k in k_range_sp_pair_expression_dist) {
  cluster_assignments_sp_pair_expression_dist[[as.character(k)]] <- cutree(row_clust_sp_pair_expression_dist, k = k)
}
silhouette_scores_sp_pair_expression_dist <- numeric(length(k_range_sp_pair_expression_dist))
for (i in seq_along(k_range_sp_pair_expression_dist)) {
  k <- k_range_sp_pair_expression_dist[i]
  sil_obj <- cluster::silhouette(cluster_assignments_sp_pair_expression_dist[[as.character(k)]], row_dist_sp_pair_expression_dist)
  silhouette_scores_sp_pair_expression_dist[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"Silhouette_scores_sp_pair_expression_dist.png"),width = ,height = ,units = ,res = 300)
plot(k_range_sp_pair_expression_dist, silhouette_scores_sp_pair_expression_dist, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
  #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_ko_pathways_sp_pair_expression_dist <- (k_range_sp_pair_expression_dist)[silhouette_scores_sp_pair_expression_dist==max(silhouette_scores_sp_pair_expression_dist)]

ko_pathways_sp_pair_expression_dist_clusters_vector <- cutree(row_clust_sp_pair_expression_dist, k = nb_clusters_ko_pathways_sp_pair_expression_dist)

  #View the row cluster assignments
print(ko_pathways_sp_pair_expression_dist_clusters_vector)

  #create pathway-to-bins correspondence vectors
v_pathway_to_sp_pair_expression_dist_clusters_or_bins <- ko_pathways_sp_pair_expression_dist_clusters_vector
names(v_pathway_to_sp_pair_expression_dist_clusters_or_bins) <- rownames(mtx_KO_pathway_specific_expression_distances)

#Chisq test between Pathway classes and pathway-specific expression distance cluster
  #contingency table
contingency_table_chisq_association_pthw_class_and_expression_dist_cluster <- table(unname(v_KO_pathway_to_most_frequent_pathway_class),unname(v_pathway_to_sp_pair_expression_dist_clusters_or_bins[names(v_KO_pathway_to_most_frequent_pathway_class)]))
contingency_table_chisq_association_pthw_class_and_expression_dist_cluster
  #pvalue
pval_chisq_association_pthw_class_and_expression_dist_cluster <- chisq.test(contingency_table_chisq_association_pthw_class_and_expression_dist_cluster)$p.value
pval_chisq_association_pthw_class_and_expression_dist_cluster
  #Cramer's V^2
squared_cramerV_association_pthw_class_and_expression_dist_cluster <- lsr::cramersV(contingency_table_chisq_association_pthw_class_and_expression_dist_cluster)^2
squared_cramerV_association_pthw_class_and_expression_dist_cluster

#Correlation between KO pathway orthogroups shared by species pairs
    #initialize the correlation matrix between ko_pathways Orthogroups correlation in species pairs and divergence time (phylogentic distance)
v_lst_all_pairs_used_for_phylo_dist <- unique(get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_pairs_of_almost_21_spcs_YPD))
v_almost_21_spcs_to_their_names_with_underscore <- gsub(pattern = " ",replacement = "_",x = rownames(mtx_pairs_of_almost_21_spcs_YPD))
names(v_almost_21_spcs_to_their_names_with_underscore) <- rownames(mtx_pairs_of_almost_21_spcs_YPD)
mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs <- matrix(NA,nrow = length(v_lst_KO_pathways_shared_by_YPD_species_pairs),ncol = length(v_lst_all_pairs_used_for_phylo_dist))
rownames(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs) <- v_lst_KO_pathways_shared_by_YPD_species_pairs
colnames(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs) <- v_lst_all_pairs_used_for_phylo_dist
    #create a list to quickly retrieve the orthogroups of the different ko pathways of a particular species
lst_YPD_almost_21_spcs_ko_pathways_orthogroups_expr_data <- list()
for (the_spc in names(v_almost_21_spcs_to_their_names_with_underscore)){
  lst_YPD_almost_21_spcs_ko_pathways_orthogroups_expr_data[[the_spc]] <-  df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY %>%
    group_by(ko_pathway,species_smpl_lbl,Orthogroup) %>%
    summarise(expression=mean(log2_1p_FPKM,na.rm=T)) %>%
    filter(species_smpl_lbl==the_spc)
}
    #calculate each correlation
for (i in 1:nrow(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs)){
  current_ko_pathway <- rownames(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs)[i]
  for (j in 1:ncol(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs)){
    #species names WITHOUT and WITH the underscore character ("_")
    current_species_pair_without_underscore_in_names <- colnames(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs)[j]
    sp1_name_without_underscore <- unlist(strsplit(current_species_pair_without_underscore_in_names,split=";"))[1]
    sp2_name_without_underscore <- unlist(strsplit(current_species_pair_without_underscore_in_names,split=";"))[2]
    sp1_name_with_underscore <- v_almost_21_spcs_to_their_names_with_underscore[sp1_name_without_underscore]
    sp2_name_with_underscore <- v_almost_21_spcs_to_their_names_with_underscore[sp2_name_without_underscore]
    #calculate current correlation as R^2
    df_current_ko_pathway_data_in_sp1 <- lst_YPD_almost_21_spcs_ko_pathways_orthogroups_expr_data[[sp1_name_without_underscore]] %>%
      filter((ko_pathway==current_ko_pathway)) %>% as.data.frame()
    if (nrow(df_current_ko_pathway_data_in_sp1)==0){
      mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs[i,j] <- NA
      next()
    }
    df_current_ko_pathway_data_in_sp2 <- lst_YPD_almost_21_spcs_ko_pathways_orthogroups_expr_data[[sp2_name_without_underscore]] %>%
      filter((ko_pathway==current_ko_pathway) & (Orthogroup %in% df_current_ko_pathway_data_in_sp1$Orthogroup)) %>% as.data.frame()
    v_common_Orthogroups_current_species_pair <- sort(intersect(df_current_ko_pathway_data_in_sp1$Orthogroup,df_current_ko_pathway_data_in_sp2$Orthogroup))
    
    rownames(df_current_ko_pathway_data_in_sp2) <- df_current_ko_pathway_data_in_sp2$Orthogroup
    df_current_ko_pathway_data_in_sp2 <- df_current_ko_pathway_data_in_sp2[v_common_Orthogroups_current_species_pair,]
    df_current_ko_pathway_data_in_sp1 <- subset(df_current_ko_pathway_data_in_sp1,Orthogroup %in% v_common_Orthogroups_current_species_pair)
    rownames(df_current_ko_pathway_data_in_sp1) <- df_current_ko_pathway_data_in_sp1$Orthogroup
    df_current_ko_pathway_data_in_sp1 <- df_current_ko_pathway_data_in_sp1[rownames(df_current_ko_pathway_data_in_sp2),]
    
    if (nrow(df_current_ko_pathway_data_in_sp2)<=2){ #if you cannot make a regression because there are not enough points to do so (at least 2)
      mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs[i,j] <- NA
      next()
    }
    
    mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs[i,j] <- (cor(df_current_ko_pathway_data_in_sp1$expression,df_current_ko_pathway_data_in_sp2$expression))^2
    
  }
  print(paste0(i," ko_pathway(s) done out of ",nrow(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs),"!"))
}

#random sampling function for power analysis (small effect size in expression distance differences)
rx_pathway <- function(n) { rlnorm(n, meanlog = 0, sdlog = 0.06)}
ry_background <- function(n) { rlnorm(n, meanlog = log(1.12), sdlog = 0.06)} #log(1.12)/0.06 = 0.2 = small Cohen's d

#compare the pathway-specific expression distances between species pairs to distances based on random sets of genes of the same size from other pathways (independent sample)
df_ko_pathways_evo_force <- data.frame(ko_pathway=v_lst_KO_pathways_shared_by_YPD_species_pairs,freq_pathway_distances_less_than_random_distances=NA,freq_pathway_distances_greater_than_random_distances=NA,
                                       freq_pathway_distances_equal_random_distances=NA,most_likely_evo_force=NA,type_evo_force=NA,pval_evo_force_call=NA,linear_rsq=NA,linear_pval=NA,
                                       has_significant_positive_parabolic_pattern=NA,min_size_pathway_in_spcs_pairs = NA, mean_size_pathway_in_spcs_pairs = NA, max_size_pathway_in_spcs_pairs = NA,
                                       sd_size_pathway_in_spcs_pairs = NA, nb_species_pairs_sharing_pathway = NA, ko_pathway_cluster = NA, power_to_detect_minimally_relevant_corr_shared_orthogroups_corr_and_divergence_time = NA,
                                       power_to_detect_minimally_relevant_diff_between_ko_pathway_dists_and_random_geneset_dists = NA, avg_slope_cor_vs_div_null_model = NA, sd_slope_cor_vs_div_null_model = NA,
                                       pval_one_sample_ttest_slope_cor_vs_div_current_pathway_greater_than_null_model = NA, stringsAsFactors = F)
the_sp_pairs <- names(mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs[the_ko_pthwy,])
for (i in 1:nrow(df_ko_pathways_evo_force)){
  the_ko_pthwy <- df_ko_pathways_evo_force$ko_pathway[i]
  
  v_current_pathway_size_in_sp_pairs <- as.vector(mtx_KO_pathway_size_in_sp_pairs[the_ko_pthwy,])
  #v_current_pathway_linreg_power_in_sp_pairs <- vapply(X = v_current_pathway_size_in_sp_pairs,FUN = function(the_smpl_size) pwr::pwr.f2.test(u = 1, v=the_smpl_size-2, f2 = 0.1, sig.level = 0.05, power = NULL)$power,FUN.VALUE = 0.0)
  df_ko_pathways_evo_force$min_size_pathway_in_spcs_pairs[i] <- min(v_current_pathway_size_in_sp_pairs,na.rm=T)
  df_ko_pathways_evo_force$mean_size_pathway_in_spcs_pairs[i] <- mean(v_current_pathway_size_in_sp_pairs,na.rm=T)
  df_ko_pathways_evo_force$max_size_pathway_in_spcs_pairs[i] <- max(v_current_pathway_size_in_sp_pairs,na.rm=T)
  df_ko_pathways_evo_force$sd_size_pathway_in_spcs_pairs[i] <- sd(v_current_pathway_size_in_sp_pairs,na.rm=T)
  df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i] <- sum((v_current_pathway_size_in_sp_pairs>2)&(!is.na(v_current_pathway_size_in_sp_pairs)),na.rm = T)
  df_ko_pathways_evo_force$ko_pathway_cluster[i] <- unname(v_pathway_to_size_clusters_or_bins[the_ko_pthwy])
  #df_ko_pathways_evo_force$min_power_orthogroups_lin_regr[i] <- min(v_current_pathway_linreg_power_in_sp_pairs,na.rm = T)
  #df_ko_pathways_evo_force$mean_power_orthogroups_lin_regr[i] <- mean(v_current_pathway_linreg_power_in_sp_pairs,na.rm = T)
  #df_ko_pathways_evo_force$median_power_orthogroups_lin_regr[i] <- median(v_current_pathway_linreg_power_in_sp_pairs,na.rm = T)
  #df_ko_pathways_evo_force$max_power_orthogroups_lin_regr[i] <- max(v_current_pathway_linreg_power_in_sp_pairs,na.rm = T)
  #df_ko_pathways_evo_force$sd_power_orthogroups_lin_regr[i] <- sd(v_current_pathway_linreg_power_in_sp_pairs,na.rm = T)
  #if the pathway does not have Orthogroups occurence in more than 2 species pair, do not make a call 
  if (df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i]<=2){ #if the sample size (number of species pairs in which the correlation between orthologs from the pathway could be assessed) is not high enough to make a regression or find a reliable one
    print(paste0(i," ko_pathway(s) done out of ",nrow(df_ko_pathways_evo_force),"!"))
    next()
  }
  df_ko_pathways_evo_force$power_to_detect_minimally_relevant_corr_shared_orthogroups_corr_and_divergence_time[i] <- pwr::pwr.f2.test(u = 1, v=df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i]-2, f2 = 0.1, sig.level = 0.05, power = NULL)$power
  df_ko_pathways_evo_force$power_to_detect_minimally_relevant_diff_between_ko_pathway_dists_and_random_geneset_dists[i] <- sum(MKpower::sim.power.wilcox.test(nx = df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i],rx = rx_pathway,rx.H0 = ry_background,
                                                                                                                                                              ny = df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i] , ry = ry_background, ry.H0 = ry_background, sig.level = 0.05,conf.int=F,alternative = "two.sided",approximate=F,iter =1000)$Exact$H1$pvalue<0.05)/1000
  
  current_pathway_df_corr_model <- data.frame(sp_pair=the_sp_pairs,sp1_name_without_underscore=vapply(X = the_sp_pairs,FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,split=";"))[1],FUN.VALUE = ""),sp2_name_without_underscore=vapply(X = the_sp_pairs,FUN = function(the_sp_pair) unlist(strsplit(the_sp_pair,split=";"))[2],FUN.VALUE = ""),x = v_phylo_dist_of_almost_21_spcs_YPD,y = mtx_cor_KO_pathway_orthogroups_expression_across_YPD_almost_21_spcs_pairs[the_ko_pthwy,],loess_prediction = NA, loess_prediction_current_change = NA, diff_x = NA, loess_preds_change_rate = NA, lbl_change_loess_pred = NA)
  current_pathway_df_corr_model <- current_pathway_df_corr_model %>% arrange(x)
  current_pathway_df_corr_model$order_sp1 <- v_order_almost_21_spcs_YPD[current_pathway_df_corr_model$sp1_name_without_underscore]
  current_pathway_df_corr_model$order_sp2 <- v_order_almost_21_spcs_YPD[current_pathway_df_corr_model$sp2_name_without_underscore]
  #if the pathway does not have Orthogroups occurence in MORE than 2 species pair, do not make a call 
  if (sum(!is.na(current_pathway_df_corr_model$y))<=2){ #if the sample size (number of species pairs in which the correlation between orthologs from the pathway could be assessed) is not high enough to make a regression or find a reliable one
    print(paste0(i," ko_pathway(s) done out of ",nrow(df_ko_pathways_evo_force),"!"))
    next()
  }

  lm_current_ko_pathways_orthogroups_correlation_across_YPD_almost_21_spcs_pairs_vs_phylo_dist <- lm(data = current_pathway_df_corr_model,formula = y~x)
  current_pathway_slope_cor_vs_div <- coef(summary(lm_current_ko_pathways_orthogroups_correlation_across_YPD_almost_21_spcs_pairs_vs_phylo_dist))[2,"Estimate"]
  current_pathway_intercept_cor_vs_div <- coef(summary(lm_current_ko_pathways_orthogroups_correlation_across_YPD_almost_21_spcs_pairs_vs_phylo_dist))["(Intercept)","Estimate"]
  current_ko_pathway_orthologs_correlation_with_divergence_residuals <- summary(lm_current_ko_pathways_orthogroups_correlation_across_YPD_almost_21_spcs_pairs_vs_phylo_dist)$resid
  current_ko_pathway_orthologs_correlation_with_divergence_boxcoxTransform_coefficient <- EnvStats::boxcox( current_ko_pathway_orthologs_correlation_with_divergence_residuals + 2*abs(min(current_ko_pathway_orthologs_correlation_with_divergence_residuals,na.rm = T)),optimize=T)$lambda
  
  #compare to random orthgroups of the same size using species pairwise distances
  v_current_ko_pathway_Orthogroups <- unique(subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,ko_pathway==the_ko_pthwy)$Orthogroup)
  count_permutations_freq_pathway_distances_equal_random_distances <- 0
  count_permutations_freq_pathway_distances_greater_than_random_distances <- 0
  count_permutations_pathway_distances_less_than_random_distances <- 0
  v_slopes_null_model_cor_vs_div_for_current_pathway <- rep(NA,nb_permutations)
  for (id_perm_random_orthogroup in 1:nb_permutations){
    #generate a set of random orthogroups of the same size as the current pathway
      #Compare the distances in that random pathway to the distances observed in the current pathway
    v_random_set_of_Orthogroups <- sample(x = setdiff(c(v_lst_accessory_Orthogroups_YPD,v_lst_core_Orthogroups_YPD),v_current_ko_pathway_Orthogroups),size = length(v_current_ko_pathway_Orthogroups) ,replace = F)
    v_null_model_orthogroups_dist_for_current_pathway <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(vegan::vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_random_set_of_Orthogroups], method = "bray")))
    if (wilcox.test(mtx_KO_pathway_specific_expression_distances[the_ko_pthwy,],v_null_model_orthogroups_dist_for_current_pathway,alternative="two.sided",paired=F)$p.value>=0.05){#is equal true
      count_permutations_freq_pathway_distances_equal_random_distances <- count_permutations_freq_pathway_distances_equal_random_distances + 1
    }else{
      if (wilcox.test(mtx_KO_pathway_specific_expression_distances[the_ko_pthwy,],v_null_model_orthogroups_dist_for_current_pathway,alternative="greater",paired=F)$p.value<0.05){#is greater true
        count_permutations_freq_pathway_distances_greater_than_random_distances <- count_permutations_freq_pathway_distances_greater_than_random_distances + 1
      }else{#is less true
        count_permutations_pathway_distances_less_than_random_distances <- count_permutations_pathway_distances_less_than_random_distances + 1
      }
    }
    
    
    #retrieve the expression correlation and phylo dist data for a random pathway (set of random orthogroups) in a random set of species pairs
    v_random_species_pairs_id <- sample(1:(choose(nrow(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD),2)),size=df_ko_pathways_evo_force$nb_species_pairs_sharing_pathway[i],replace = F)
    v_expr_cor_between_orthogroups_across_sp_pairs_in_null_model_pathway <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,sample(1:ncol(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD),size = length(v_current_ko_pathway_Orthogroups),replace = F)]))^2))[v_random_species_pairs_id]
    while (any(is.na(v_expr_cor_between_orthogroups_across_sp_pairs_in_null_model_pathway))){
      v_expr_cor_between_orthogroups_across_sp_pairs_in_null_model_pathway <- get_v_uppr_tri_no_diag_from_mtx(the_mtx = as.matrix(cor(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,sample(1:ncol(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD),size = length(v_current_ko_pathway_Orthogroups),replace = F)]))^2))[v_random_species_pairs_id]
    }
    v_phylo_dist_between_sp_pairs_in_null_model_pathway <- v_phylo_dist_of_almost_21_spcs_YPD[v_random_species_pairs_id]
    v_slopes_null_model_cor_vs_div_for_current_pathway[id_perm_random_orthogroup] <- coef(summary(lm(v_expr_cor_between_orthogroups_across_sp_pairs_in_null_model_pathway~v_phylo_dist_between_sp_pairs_in_null_model_pathway)))[2,"Estimate"]
    
  }
  df_ko_pathways_evo_force$freq_pathway_distances_equal_random_distances[i] <- count_permutations_freq_pathway_distances_equal_random_distances/nb_permutations
  df_ko_pathways_evo_force$freq_pathway_distances_greater_than_random_distances[i] <- count_permutations_freq_pathway_distances_greater_than_random_distances/nb_permutations
  df_ko_pathways_evo_force$freq_pathway_distances_less_than_random_distances[i] <- count_permutations_pathway_distances_less_than_random_distances/nb_permutations
  df_ko_pathways_evo_force$avg_slope_cor_vs_div_null_model[i] <- mean(v_slopes_null_model_cor_vs_div_for_current_pathway,na.rm=T)
  df_ko_pathways_evo_force$sd_slope_cor_vs_div_null_model[i] <- sd(v_slopes_null_model_cor_vs_div_for_current_pathway,na.rm=T)
  df_ko_pathways_evo_force$pval_one_sample_ttest_slope_cor_vs_div_current_pathway_greater_than_null_model[i] <- t.test(v_slopes_null_model_cor_vs_div_for_current_pathway,mu=current_pathway_slope_cor_vs_div,alternative = "greater")$p.value
  
  #draw plot current pathway orthogroups correlation vs divergence time across species pairs sharing these orthogroups + Null Model
  ggplot(data = current_pathway_df_corr_model, aes(x = x, y = y)) +
    geom_point(aes(color = order_sp1, fill = order_sp2),shape=21,size=3,stroke=2) +
    geom_smooth(method = "loess", se = TRUE, color = "red", size = 1) + 
    geom_abline(slope = mean(v_slopes_null_model_cor_vs_div_for_current_pathway,na.rm=T), intercept = current_pathway_intercept_cor_vs_div, linetype = "dotted", color = "black",linewidth = 1.5)+ #this is the null model (slope is determine by random set of orthogroups)
    xlab("Divergence time (phylogenetic distance between a species pair)")+
    ylab("Pathway Orthogroups expression correlation in the species pair (R2)")+
    labs(title = paste0(the_ko_pthwy,":",v_KO_pathway_to_pathway_title[the_ko_pthwy]),fill = "Order")+
    scale_fill_manual(values=v_clade_color) + 
    scale_colour_manual(values=v_clade_color)+
    guides(color = "none")
  ggsave(filename = paste0(the_ko_pthwy,"_Orthogroups_correlation_in_sp_pairs_VS_divergence_time.svg"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,device = svg)
  ggsave(filename = paste0(the_ko_pthwy,"_Orthogroups_correlation_in_sp_pairs_VS_divergence_time.png"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,dpi=300)
  
  #glvma test for linear relationship p-value
  #pval_lm_current_ko_pathway_corr <- gvlma::gvlma(lm_current_ko_pathways_orthogroups_correlation_across_YPD_almost_21_spcs_pairs_vs_phylo_dist)$GlobalTest$GlobalStat4$pvalue[1,1]
  #df_ko_pathways_evo_force$lm_diagnostic_pval[i] <- pval_lm_current_ko_pathway_corr  #heuristic_cond_lm <- (((current_ko_pathway_orthologs_correlation_with_divergence_boxcoxTransform_coefficient < 1) & (abs(current_ko_pathway_orthologs_correlation_with_divergence_boxcoxTransform_coefficient-1)<0.25 ))||((current_ko_pathway_orthologs_correlation_with_divergence_boxcoxTransform_coefficient >= 1) & (abs(current_ko_pathway_orthologs_correlation_with_divergence_boxcoxTransform_coefficient-1)<0.5 )))&(max(abs(range(current_pathway_df_corr_model$y[1:round(nrow(current_pathway_df_corr_model)/10,0)])-range(current_pathway_df_corr_model$y[(round(nrow(current_pathway_df_corr_model)/1.111,0)+1):nrow(current_pathway_df_corr_model)]))) >= 0.1)
  #if ((pval_lm_current_ko_pathway_corr>=0.05)&(heuristic_cond_lm)) {}else{}
  
  #determine if the curve has a segment with a significant positive parabolic pattern that represents positive selection
  current_pathway_df_corr_model <- subset(current_pathway_df_corr_model,!is.na(current_pathway_df_corr_model$y))
  current_pathway_df_corr_model$loess_prediction <- unname(predict(loess(data = current_pathway_df_corr_model,y~x), se = TRUE)$fit)
  current_pathway_df_corr_model$loess_prediction_current_change <- c(diff(current_pathway_df_corr_model$loess_prediction),NA)
  current_pathway_df_corr_model$diff_x <- c(diff(current_pathway_df_corr_model$x),NA)
  current_pathway_df_corr_model$loess_preds_change_rate <- current_pathway_df_corr_model$loess_prediction_current_change/current_pathway_df_corr_model$diff_x
  current_pathway_df_corr_model$lbl_change_loess_pred <- ifelse(test = !is.na(current_pathway_df_corr_model$loess_preds_change_rate),yes = ifelse(test = current_pathway_df_corr_model$loess_preds_change_rate>=0.05,yes = "increasing",no = ifelse(test = current_pathway_df_corr_model$loess_preds_change_rate<=-0.05,yes = "decreasing",no = "stagnating")),no = current_pathway_df_corr_model$loess_preds_change_rate)
  df_current_pathway_loess_lbl_change_streaks <- data.frame(len_streak=rle(current_pathway_df_corr_model$lbl_change_loess_pred)$lengths,associated_change=rle(current_pathway_df_corr_model$lbl_change_loess_pred)$values,stringsAsFactors = F)
  df_current_pathway_loess_lbl_change_streaks$cumul_pt_index <- cumsum(df_current_pathway_loess_lbl_change_streaks$len_streak)
  #find considerable decreasing-stagnating-increasing or decreasing-increasing patterns (positive parabolic pattern)
    #initialize boolean
  df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- F
    #determine whether or not there is an impactful decreasing-stagnating-increasing pattern
  longVec <- df_current_pathway_loess_lbl_change_streaks$associated_change
  shortVec <- c("decreasing","stagnating","increasing")
  
      #Find potential starting points where the first element of shortVec matches
  idx <- which(longVec == shortVec[1])
  
      #Filter these potential starting points to find actual matches of the entire shortVec
  start_index_dec_stag_inc <- idx[sapply(idx, function(i) {all(longVec[i:(i + length(shortVec) - 1)] == shortVec)})][1]
  if (any(is.na(start_index_dec_stag_inc))){ #if there are no decreasing-stagnating-increasing patterns 
    df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- F
    #next()
  }else{ #there is at least one decreasing-stagnating-increasing 
    #make sure that there is at least one positive parabolic pattern that is significant in any of the segments it has been identified
    for (current_start_index_dec_stag_inc in start_index_dec_stag_inc){
      current_pos_parab_pattern_pt_index_start <- df_current_pathway_loess_lbl_change_streaks$cumul_pt_index[current_start_index_dec_stag_inc] - df_current_pathway_loess_lbl_change_streaks$len_streak[current_start_index_dec_stag_inc] + 1
      current_pos_parab_pattern_pt_index_end <- df_current_pathway_loess_lbl_change_streaks$cumul_pt_index[current_start_index_dec_stag_inc+2]
      #fit second-order model
      current_parab_fit <- lm(y ~ x + I(x^2), data = current_pathway_df_corr_model[current_pos_parab_pattern_pt_index_start:current_pos_parab_pattern_pt_index_end,])
        #get second-order x^2 coefficient "a"
      a_x2_coef_current_parab_model <- unname(coef(current_parab_fit)["I(x^2)"])
        #get model p-value
      fstats_current_parab_model <- summary(current_parab_fit)$fstatistic
      pval_current_parab_model <- pf(fstats_current_parab_model[1],fstats_current_parab_model[2],fstats_current_parab_model[3],lower.tail = F)
      if ((is.nan(unname(pval_current_parab_model)))|(nrow(current_pathway_df_corr_model[current_pos_parab_pattern_pt_index_start:current_pos_parab_pattern_pt_index_end,])<=2)){ #if the sample size is too small for a reliable regression fit
        next()
      }
      if ((a_x2_coef_current_parab_model>0)&(pval_current_parab_model<0.05)){
        df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- T
        df_ko_pathways_evo_force$most_likely_evo_force[i] <- "Positive selection"
        df_ko_pathways_evo_force$type_evo_force[i] <- "Convergent"
        df_ko_pathways_evo_force$pval_evo_force_call[i] <- 1 - df_ko_pathways_evo_force$freq_pathway_distances_less_than_random_distances[i]
        df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- T
        
        #draw plot current pathway orthogroups distance vs distribution of the distance for a set of random orthogroups of the same size (pathway Null Model)
        ggplot() +
          geom_density(aes(mtx_KO_pathway_specific_expression_distances[the_ko_pthwy,],col="Observed",fill="Observed"),alpha=0.3) +
          geom_density(aes(v_null_model_orthogroups_dist_for_current_pathway,col="Null Model",fill="Null Model"),alpha=0.3) +
          xlab("Pathway-specific expression distance between species pairs")+
          ylab("Density")+
          labs(title = paste0(the_ko_pthwy,":",v_KO_pathway_to_pathway_title[the_ko_pthwy],"\n(Evolutionary force call: ",df_ko_pathways_evo_force$most_likely_evo_force[i],"; p-value: ",df_ko_pathways_evo_force$pval_evo_force_call[i],")"))+
          scale_fill_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name = "Pathway-specific\nSpecies distances") + 
          scale_colour_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name="")+
          guides(color = "none")
        ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.svg"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,device = svg)
        ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.png"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,dpi=300)
        
        break()
      }
    }
  }
  
  #at this stage, determine whether or not there is an impactful decreasing-increasing pattern IF positive selection has not been detected yet
  if (!df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i]){
    longVec <- df_current_pathway_loess_lbl_change_streaks$associated_change
    shortVec <- c("decreasing","increasing")
    
    #Find potential starting points where the first element of shortVec matches
    idx <- which(longVec == shortVec[1])
    
    #Filter these potential starting points to find actual matches of the entire shortVec
    start_index_dec_inc <- idx[sapply(idx, function(i) {all(longVec[i:(i + length(shortVec) - 1)] == shortVec)})][1]
    if (any(is.na(start_index_dec_inc))){ #if there are no decreasing-increasing patterns 
      df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- F
      #next()
    }else{ #there is at least one decreasing-increasing 
      #make sure that there is at least one positive parabolic pattern that is significant in any of the segments it has been identified
      for (current_start_index_dec_inc in start_index_dec_inc){
        current_pos_parab_pattern_pt_index_start <- df_current_pathway_loess_lbl_change_streaks$cumul_pt_index[current_start_index_dec_inc] - df_current_pathway_loess_lbl_change_streaks$len_streak[current_start_index_dec_inc] + 1
        current_pos_parab_pattern_pt_index_end <- df_current_pathway_loess_lbl_change_streaks$cumul_pt_index[current_start_index_dec_inc+1]
        #fit second-order model
        current_parab_fit <- lm(y ~ x + I(x^2), data = current_pathway_df_corr_model[current_pos_parab_pattern_pt_index_start:current_pos_parab_pattern_pt_index_end,])
        #get second-order x^2 coefficient "a"
        a_x2_coef_current_parab_model <- unname(coef(current_parab_fit)["I(x^2)"])
        #get model p-value
        fstats_current_parab_model <- summary(current_parab_fit)$fstatistic
        pval_current_parab_model <- pf(fstats_current_parab_model[1],fstats_current_parab_model[2],fstats_current_parab_model[3],lower.tail = F)
        if ((is.nan(unname(pval_current_parab_model)))|(nrow(current_pathway_df_corr_model[current_pos_parab_pattern_pt_index_start:current_pos_parab_pattern_pt_index_end,])<=2)){ #if the sample size is too small for a reliable regression fit
          next()
        }
        if ((a_x2_coef_current_parab_model>0)&(pval_current_parab_model<0.05)){
          df_ko_pathways_evo_force$has_significant_positive_parabolic_pattern[i] <- T
          df_ko_pathways_evo_force$most_likely_evo_force[i] <- "Positive selection"
          df_ko_pathways_evo_force$type_evo_force[i] <- "Convergent"
          df_ko_pathways_evo_force$pval_evo_force_call[i] <- 1 - df_ko_pathways_evo_force$freq_pathway_distances_less_than_random_distances[i]
          
          #draw plot current pathway orthogroups distance vs distribution of the distance for a set of random orthogroups of the same size (pathway Null Model)
          ggplot() +
            geom_density(aes(mtx_KO_pathway_specific_expression_distances[the_ko_pthwy,],col="Observed",fill="Observed"),alpha=0.3) +
            geom_density(aes(v_null_model_orthogroups_dist_for_current_pathway,col="Null Model",fill="Null Model"),alpha=0.3) +
            xlab("Pathway-specific expression distance between species pairs")+
            ylab("Density")+
            labs(title = paste0(the_ko_pthwy,":",v_KO_pathway_to_pathway_title[the_ko_pthwy],"\n(Evolutionary force call: ",df_ko_pathways_evo_force$most_likely_evo_force[i],"; p-value: ",df_ko_pathways_evo_force$pval_evo_force_call[i],")"))+
            scale_fill_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name = "Pathway-specific\nSpecies distances") + 
            scale_colour_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name="")+
            guides(color = "none")
          ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.svg"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,device = svg)
          ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.png"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,dpi=300)
          
          break()
        }
      }
    }
  }
  #if positive selection has been detected, go to the next iteration
  if ((!is.na(df_ko_pathways_evo_force$most_likely_evo_force[i]))&(df_ko_pathways_evo_force$most_likely_evo_force[i]=="Positive selection")){
    print(paste0(i," ko_pathway(s) done out of ",nrow(df_ko_pathways_evo_force),"!"))
    next()
  }
  #if the code below did not go to the next iteration, then we can discard positive selection as being the predominant evolutionary force for the current pathway
  current_max_freq <- max(c(df_ko_pathways_evo_force$freq_pathway_distances_equal_random_distances[i],df_ko_pathways_evo_force$freq_pathway_distances_greater_than_random_distances[i],df_ko_pathways_evo_force$freq_pathway_distances_less_than_random_distances[i]),na.rm = T)
  df_ko_pathways_evo_force$pval_evo_force_call[i] <- 1 - current_max_freq
  df_ko_pathways_evo_force$most_likely_evo_force[i] <- c("Neutral","Diversifying selection","Negative selection")[which(c(df_ko_pathways_evo_force$freq_pathway_distances_equal_random_distances[i],df_ko_pathways_evo_force$freq_pathway_distances_greater_than_random_distances[i],df_ko_pathways_evo_force$freq_pathway_distances_less_than_random_distances[i])==current_max_freq)]
  df_ko_pathways_evo_force$type_evo_force[i] <- ifelse(test = df_ko_pathways_evo_force$most_likely_evo_force[i]%in%c("Neutral","Diversifying selection"),yes = "Divergent",no = "Convergent")
  
  #draw plot current pathway orthogroups distance vs distribution of the distance for a set of random orthogroups of the same size (pathway Null Model)
  ggplot() +
    geom_density(aes(mtx_KO_pathway_specific_expression_distances[the_ko_pthwy,],col="Observed",fill="Observed"),alpha=0.3) +
    geom_density(aes(v_null_model_orthogroups_dist_for_current_pathway,col="Null Model",fill="Null Model"),alpha=0.3) +
    xlab("Pathway-specific expression distance between species pairs")+
    ylab("Density")+
    labs(title = paste0(the_ko_pthwy,":",v_KO_pathway_to_pathway_title[the_ko_pthwy],"\n(Evolutionary force call: ",df_ko_pathways_evo_force$most_likely_evo_force[i],"; p-value: ",df_ko_pathways_evo_force$pval_evo_force_call[i],")"))+
    scale_fill_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name = "Pathway-specific\nSpecies distances") + 
    scale_colour_manual(values=c("Observed"="dodgerblue","Null Model"="grey70"), name="")+
    guides(color = "none")
  ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.svg"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,device = svg)
  ggsave(filename = paste0(the_ko_pthwy,"_Distribution_pathwaySpecific_expr_distance_vs_NullModel.png"), path=paste0(output_workspace,"ko_pathway_plots/"), width = 15/2.54, height = 15/2.54,dpi=300)
       
  print(paste0(i," ko_pathway(s) done out of ",nrow(df_ko_pathways_evo_force),"!"))
  
}
df_ko_pathways_evo_force$pval_evo_force_call <- p.adjust(df_ko_pathways_evo_force$pval_evo_force_call,"fdr")
df_ko_pathways_evo_force$pathway_title <- v_KO_pathway_to_pathway_title[df_ko_pathways_evo_force$ko_pathway]
df_ko_pathways_evo_force$most_frequent_pathway_class <- v_KO_pathway_to_most_frequent_pathway_class[df_ko_pathways_evo_force$ko_pathway]
df_ko_pathways_evo_force$most_frequent_pathway_group <- v_KO_pathway_to_most_frequent_pathway_group[df_ko_pathways_evo_force$ko_pathway]

#create 2 power bins (one that determines if power >= 0.8 for the pathway linear regression; the other determines if power >= 0.8 for the pathway-based pairwise distance comparisons)
df_ko_pathways_evo_force$is_power_enough_for_pathway_linreg <- df_ko_pathways_evo_force$power_to_detect_minimally_relevant_corr_shared_orthogroups_corr_and_divergence_time >= 0.8
df_ko_pathways_evo_force$is_power_enough_for_pathway_linreg[is.na(df_ko_pathways_evo_force$is_power_enough_for_pathway_linreg)] <- F
df_ko_pathways_evo_force$is_power_enough_for_sp_pairs_pathway_based_dist_comparison <- df_ko_pathways_evo_force$power_to_detect_minimally_relevant_diff_between_ko_pathway_dists_and_random_geneset_dists >= 0.8
df_ko_pathways_evo_force$is_power_enough_for_sp_pairs_pathway_based_dist_comparison[is.na(df_ko_pathways_evo_force$is_power_enough_for_sp_pairs_pathway_based_dist_comparison)] <- F

#summarize by species pair prevalence clusters of orthgroups
df_pathway_bins_evo_force <- df_ko_pathways_evo_force %>% 
  filter(!is.na(most_likely_evo_force)) %>%
  group_by(ko_pathway_cluster) %>%
  summarise(avg_nb_species_pairs_in_which_the_pathway_is_shared=mean(unname(rowSums((mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),]>2)&(!is.na(mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),])))),na.rm=T),
            avg_power_based_on_nb_species_pairs=mean(vapply(X = unname(rowSums((mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),]>2)&(!is.na(mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),])))),FUN = function(the_smpl_size) pwr::pwr.f2.test(u = 1, v=the_smpl_size-2, f2 = 0.1, sig.level = 0.05, power = NULL)$power,FUN.VALUE = 0.0),na.rm=T),
            sd_power_based_on_nb_species_pairs=sd(vapply(X = unname(rowSums((mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),]>2)&(!is.na(mtx_KO_pathway_size_in_sp_pairs[unique(ko_pathway),])))),FUN = function(the_smpl_size) pwr::pwr.f2.test(u = 1, v=the_smpl_size-2, f2 = 0.1, sig.level = 0.05, power = NULL)$power,FUN.VALUE = 0.0),na.rm=T),
            freq_pos_sel=sum((most_likely_evo_force=="Positive selection")&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_neg_sel=sum((most_likely_evo_force=="Negative selection")&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_conv_forces=sum((most_likely_evo_force%in%c("Positive selection","Negative selection"))&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_neutral=sum((most_likely_evo_force=="Neutral")&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_div_sel=sum((most_likely_evo_force=="Diversifying selection")&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_div_forces=sum((most_likely_evo_force%in%c("Neutral","Diversifying selection"))&(!is.na(most_likely_evo_force)))/length(ko_pathway),
            freq_signif_pos_sel=sum((most_likely_evo_force=="Positive selection")&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_signif_neg_sel=sum((most_likely_evo_force=="Negative selection")&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_signif_conv_forces=sum((most_likely_evo_force%in%c("Positive selection","Negative selection"))&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_signif_neutral=sum((most_likely_evo_force=="Neutral")&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_signif_div_sel=sum((most_likely_evo_force=="Diversifying selection")&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_signif_div_forces=sum((most_likely_evo_force%in%c("Neutral","Diversifying selection"))&(!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call<0.05))/length(ko_pathway),
            freq_NS_evo_force_calls=sum((!is.na(most_likely_evo_force))&(!is.na(pval_evo_force_call))&(pval_evo_force_call>=0.05))/length(ko_pathway),
            freq_not_called_pathways=sum((is.na(most_likely_evo_force))|(is.na(pval_evo_force_call)))/length(ko_pathway)) %>%
  arrange(avg_nb_species_pairs_in_which_the_pathway_is_shared)
df_pathway_bins_evo_force$new_bins_id_based_on_increasing_nb_sp_pairs <- paste0("Cluster",1:nrow(df_pathway_bins_evo_force),"\n(mean_n = ",round(df_pathway_bins_evo_force$avg_nb_species_pairs_in_which_the_pathway_is_shared,2),";\npower = ",round(df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs,2),"\u00B1",round(df_pathway_bins_evo_force$sd_power_based_on_nb_species_pairs,2),")")

  #plot the regression between avg_power of a cluster and the frequency of significant calls
    #create column with the freq of significant calls
df_pathway_bins_evo_force$freq_signif_calls <- 1-df_pathway_bins_evo_force$freq_NS_evo_force_calls
df_pathway_bins_evo_force$freq_signif_calls_p1pct <- df_pathway_bins_evo_force$freq_signif_calls + 0.01
df_pathway_bins_evo_force$log10p1pct_freq_signif_calls <- log10(df_pathway_bins_evo_force$freq_signif_calls_p1pct)
    #create exponential model with the linear model log10(y+0.01)~x (0.01 is a little frequency of 1% that act as a constant to avoid log10(0), which is NaN)
test_model_exp_signifcall_vs_avg_power <- lm(df_pathway_bins_evo_force$log10p1pct_freq_signif_calls~df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs)
    #Create predicted values vector from model (exp(predicted_log10p1pct_y))
predictions_test_model_exp_signifcall_vs_avg_power <- 10^(predict(
  test_model_exp_signifcall_vs_avg_power, 
  newdata = df_pathway_bins_evo_force[,c("avg_power_based_on_nb_species_pairs","freq_signif_calls_p1pct")],
  interval="prediction",
  level = 0.95
))
    #Plot the original data and the regression line
      #exponential fit on raw data
svg(filename = paste0(output_workspace,"RAW_Freq_signif_evo_calls_vs_avg_power_in_ko_pathway_WITH_EXPO_FIT.svg"),width = 17.4/2.54,height = 13.2/2.54)
plot(df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs, df_pathway_bins_evo_force$freq_signif_calls, main="Exponential Regression", xlab="Average power of linear regression for\northogroups of the KOG pathway", 
     ylab="Frequency of significant evolutionary calls", pch=19)
lines(df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs, predictions_test_model_exp_signifcall_vs_avg_power[,1]-0.01, col="red", lty=2)
lines(df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs, predictions_test_model_exp_signifcall_vs_avg_power[,2]-0.01, col="blue", lty=2)
lines(df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs, predictions_test_model_exp_signifcall_vs_avg_power[,3]-0.01, col="blue", lty=2)
            #legend("topright", legend="Exponential Regression", col="red", lwd=2)
dev.off()
      #linear model on log-transformed data
df_pathway_bins_evo_force$x <- df_pathway_bins_evo_force$avg_power_based_on_nb_species_pairs
df_pathway_bins_evo_force$y <- df_pathway_bins_evo_force$log10p1pct_freq_signif_calls
lbl_eqn_rsq_freq_signif_evo_calls_vs_avg_power_in_ko_pathway <- lm_eqn(df = df_pathway_bins_evo_force)
ggplot(data = df_pathway_bins_evo_force) +
  geom_point(mapping = aes(x = avg_power_based_on_nb_species_pairs, y = log10p1pct_freq_signif_calls)) + 
  geom_smooth(mapping = aes(x = avg_power_based_on_nb_species_pairs, y = log10p1pct_freq_signif_calls),method = "lm", se=FALSE, color="black",lty=2) +
  geom_text(x = 0.25, y = -0.8, label = lbl_eqn_rsq_freq_signif_evo_calls_vs_avg_power_in_ko_pathway, parse = TRUE, cex = 6) + 
  xlab("Average power of linear regression for\northogroups of the KOG pathway") + 
  ylab("log10(Frequency of significant evolutionary calls + 0.01)")+
  theme(legend.key.size = unit(1, "cm")) #+  scale_y_continuous(breaks = seq(0,1,0.1),limits = c(0,1))
ggsave(filename = "LOG-TRANSFORMED_Freq_signif_evo_calls_vs_avg_power_in_ko_pathway_WITH_EXPO_FIT.svg", path=output_workspace, width = 17.4/2.54, height = 13.2/2.54,device = svg)

  #plot the frequency of significant calls for different clusters (using new_bins_id_based_on_increasing_nb_sp_pairs)
ggplot(data = df_pathway_bins_evo_force) +
  geom_col(mapping = aes(x = new_bins_id_based_on_increasing_nb_sp_pairs, y = freq_signif_calls),fill="cornflowerblue") + 
  xlab("Pathway presence/absence cluster (average number of species sharing it)") + 
  ylab("Frequency of significant evolutionary calls") +
  theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=14))
ggsave(filename = "Freq_signif_evo_calls_vs_orthogroup_CLUSTER_based_on_prevalence_in_sp_pairs.svg", path=output_workspace, width = 30/2.54, height = 25/2.54,device = svg)

#    #create a matrix of the variance of ko pathways orthgroups within each species
#mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs <- matrix(NA,nrow = length(v_lst_KO_pathways_shared_by_YPD_species_pairs),ncol = length(v_almost_21_spcs_to_their_names_with_underscore))
#rownames(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs) <- v_lst_KO_pathways_shared_by_YPD_species_pairs
#colnames(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs) <- names(v_almost_21_spcs_to_their_names_with_underscore)
#    #calculate each correlation
#for (i in 1:nrow(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs)){
#  current_ko_pathway <- rownames(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs)[i]
#  for (j in 1:ncol(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs)){
#    #species names WITHOUT and WITH the underscore character ("_")
#    current_species_without_underscore_in_names <- colnames(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs)[j]
#    current_species_name_with_underscore <- v_almost_21_spcs_to_their_names_with_underscore[current_species_without_underscore_in_names]
#    #calculate current correlation as R^2
#    df_current_ko_pathway_data_in_current_species <- lst_YPD_almost_21_spcs_ko_pathways_orthogroups_expr_data[[current_species_without_underscore_in_names]] %>%
#      filter((ko_pathway==current_ko_pathway)) %>% as.data.frame()
#    if (nrow(df_current_ko_pathway_data_in_current_species)==0){
#      mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs[i,j] <- NA
#      next()
#    }
#    
#    
#    mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs[i,j] <- var(df_current_ko_pathway_data_in_current_species$expression,na.rm = T)
#    
#  }
#  print(paste0(i," ko_pathway(s) done out of ",nrow(mtx_KO_pathway_orthogroups_expression_variance_across_YPD_almost_21_spcs),"!"))
#}

#test whether or not a pathway class is enriched among positively selected pathways 
v_unique_ko_pathway_classes <- unique(unname(v_KO_pathway_to_most_frequent_pathway_class))
v_pvals_pathway_class_enriched_in_positively_selected_ko_pathways <- rep(NA, length(v_unique_ko_pathway_classes))
names(v_pvals_pathway_class_enriched_in_positively_selected_ko_pathways) <- v_unique_ko_pathway_classes
for (id_class in 1:length(v_unique_ko_pathway_classes)){
  current_ko_pathway_class <- v_unique_ko_pathway_classes[id_class]
  v_current_class_ko_pathways <- rownames(mtx_pathway_to_pthw_class_association)[mtx_pathway_to_pthw_class_association[,current_ko_pathway_class]>0]
  mtx_cont_current_class <- matrix(NA,2,2)
  rownames(mtx_cont_current_class) <- c(paste0("Not ",current_ko_pathway_class),current_ko_pathway_class)
  colnames(mtx_cont_current_class) <- c("Other evolutionary forces","Positive selection")
  mtx_cont_current_class[1,1] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Positive selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[2,1] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Positive selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[1,2] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Positive selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[2,2] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Positive selection"))$ko_pathway,v_current_class_ko_pathways))
  v_pvals_pathway_class_enriched_in_positively_selected_ko_pathways[id_class] <- fisher.test(mtx_cont_current_class,alternative = "greater")$p.value
}
v_adj_pvals_fdr_pathway_class_enriched_in_positively_selected_ko_pathways <- p.adjust(v_pvals_pathway_class_enriched_in_positively_selected_ko_pathways,"fdr")
v_adj_pvals_fdr_pathway_class_enriched_in_positively_selected_ko_pathways


#test whether or not a pathway group is enriched among positively selected pathways 
v_unique_ko_pathway_groups <- unique(unname(v_KO_pathway_to_most_frequent_pathway_group))
v_pvals_pathway_group_enriched_in_positively_selected_ko_pathways <- rep(NA, length(v_unique_ko_pathway_groups))
names(v_pvals_pathway_group_enriched_in_positively_selected_ko_pathways) <- v_unique_ko_pathway_groups
for (id_group in 1:length(v_unique_ko_pathway_groups)){
  current_ko_pathway_group <- v_unique_ko_pathway_groups[id_group]
  v_current_group_ko_pathways <- rownames(mtx_pathway_to_pthw_group_association)[mtx_pathway_to_pthw_group_association[,current_ko_pathway_group]>0]
  mtx_cont_current_group <- matrix(NA,2,2)
  rownames(mtx_cont_current_group) <- c(paste0("Not ",current_ko_pathway_group),current_ko_pathway_group)
  colnames(mtx_cont_current_group) <- c("Other evolutionary forces","Positive selection")
  mtx_cont_current_group[1,1] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Positive selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[2,1] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Positive selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[1,2] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Positive selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[2,2] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Positive selection"))$ko_pathway,v_current_group_ko_pathways))
  v_pvals_pathway_group_enriched_in_positively_selected_ko_pathways[id_group] <- fisher.test(mtx_cont_current_group,alternative = "greater")$p.value
}
v_adj_pvals_fdr_pathway_group_enriched_in_positively_selected_ko_pathways <- p.adjust(v_pvals_pathway_group_enriched_in_positively_selected_ko_pathways,"fdr")
v_adj_pvals_fdr_pathway_group_enriched_in_positively_selected_ko_pathways


#test whether or not a pathway class is enriched among negatively selected pathways 
v_unique_ko_pathway_classes <- unique(unname(v_KO_pathway_to_most_frequent_pathway_class))
v_pvals_pathway_class_enriched_in_negatively_selected_ko_pathways <- rep(NA, length(v_unique_ko_pathway_classes))
names(v_pvals_pathway_class_enriched_in_negatively_selected_ko_pathways) <- v_unique_ko_pathway_classes
for (id_class in 1:length(v_unique_ko_pathway_classes)){
  current_ko_pathway_class <- "Genetic Information Processing"#v_unique_ko_pathway_classes[id_class]
  v_current_class_ko_pathways <- rownames(mtx_pathway_to_pthw_class_association)[mtx_pathway_to_pthw_class_association[,current_ko_pathway_class]>0]
  mtx_cont_current_class <- matrix(NA,2,2)
  rownames(mtx_cont_current_class) <- c(paste0("Not ",current_ko_pathway_class),current_ko_pathway_class)
  colnames(mtx_cont_current_class) <- c("Other evolutionary forces","Negative selection")
  mtx_cont_current_class[1,1] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Negative selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[2,1] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Negative selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[1,2] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Negative selection"))$ko_pathway,v_current_class_ko_pathways))
  mtx_cont_current_class[2,2] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Negative selection"))$ko_pathway,v_current_class_ko_pathways))
  v_pvals_pathway_class_enriched_in_negatively_selected_ko_pathways[id_class] <- fisher.test(mtx_cont_current_class,alternative = "greater")$p.value
}
v_adj_pvals_fdr_pathway_class_enriched_in_negatively_selected_ko_pathways <- p.adjust(v_pvals_pathway_class_enriched_in_negatively_selected_ko_pathways,"fdr")
v_adj_pvals_fdr_pathway_class_enriched_in_negatively_selected_ko_pathways


#test whether or not a pathway group is enriched among negatively selected pathways 
v_unique_ko_pathway_groups <- unique(unname(v_KO_pathway_to_most_frequent_pathway_group))
v_pvals_pathway_group_enriched_in_negatively_selected_ko_pathways <- rep(NA, length(v_unique_ko_pathway_groups))
names(v_pvals_pathway_group_enriched_in_negatively_selected_ko_pathways) <- v_unique_ko_pathway_groups
for (id_group in 1:length(v_unique_ko_pathway_groups)){
  current_ko_pathway_group <- v_unique_ko_pathway_groups[id_group]
  v_current_group_ko_pathways <- rownames(mtx_pathway_to_pthw_group_association)[mtx_pathway_to_pthw_group_association[,current_ko_pathway_group]>0]
  mtx_cont_current_group <- matrix(NA,2,2)
  rownames(mtx_cont_current_group) <- c(paste0("Not ",current_ko_pathway_group),current_ko_pathway_group)
  colnames(mtx_cont_current_group) <- c("Other evolutionary forces","Negative selection")
  mtx_cont_current_group[1,1] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Negative selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[2,1] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force=="Negative selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[1,2] <- length(intersect(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Negative selection"))$ko_pathway,v_current_group_ko_pathways))
  mtx_cont_current_group[2,2] <- length(setdiff(subset(df_ko_pathways_evo_force,(!is.na(most_likely_evo_force))&(most_likely_evo_force!="Negative selection"))$ko_pathway,v_current_group_ko_pathways))
  v_pvals_pathway_group_enriched_in_negatively_selected_ko_pathways[id_group] <- fisher.test(mtx_cont_current_group,alternative = "greater")$p.value
}
v_adj_pvals_fdr_pathway_group_enriched_in_negatively_selected_ko_pathways <- p.adjust(v_pvals_pathway_group_enriched_in_negatively_selected_ko_pathways,"fdr")
v_adj_pvals_fdr_pathway_group_enriched_in_negatively_selected_ko_pathways

#correct for multiple testing
df_ko_pathways_evo_force$most_likely_and_signif_evo_force <- ifelse(test = !is.na(df_ko_pathways_evo_force$pval_evo_force_call)&(df_ko_pathways_evo_force$pval_evo_force_call<0.05),yes = df_ko_pathways_evo_force$most_likely_evo_force,no = NA)
df_ko_pathways_evo_force$signif_type_evo_force <- ifelse(test = !is.na(df_ko_pathways_evo_force$pval_evo_force_call)&(df_ko_pathways_evo_force$pval_evo_force_call<0.05),yes = df_ko_pathways_evo_force$type_evo_force,no = NA)
df_ko_pathways_evo_force$is_evo_force_call_significant <- df_ko_pathways_evo_force$pval_evo_force_call < 0.05

#Illustrate the frequency of evolutionary convergence across power bins using a ggplot pie chart
  #Make grid of piecharts
plot_data <- df_ko_pathways_evo_force %>%
  count(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison, signif_type_evo_force) %>%
  group_by(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison) %>%
  mutate(percent = n / sum(n))

ggplot(plot_data, aes(x = "", y = percent, fill = signif_type_evo_force)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(is_power_enough_for_pathway_linreg ~ is_power_enough_for_sp_pairs_pathway_based_dist_comparison) +
  xlab("Is the power of the Mann-Whitney test >0.8?")+
  ylab("Is the power of the linear regression >0.8?")+
  labs(title = "", fill = "Type of evolutionary force acting at\nthe transcriptomic level for the KEGG pathways")
ggsave(filename = "Piechart_evo_force_type_vs_power_bins_KEGG_pathways.svg", path=output_workspace, width = 17.4/2.54, height = 13.2/2.54,device = svg)
ggsave(filename = "Piechart_evo_force_type_vs_power_bins_KEGG_pathways.png", path=output_workspace, width = 17.4/2.54, height = 13.2/2.54,device = png,dpi = 300)
 #make it as a barplot
barplot_data <- df_ko_pathways_evo_force %>%
  count(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison, is_evo_force_call_significant, type_evo_force) %>%
  group_by(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison,is_evo_force_call_significant) %>%
  mutate(percent = n / sum(n))

barplot_data$is_evo_force_call_significant <- ifelse(test = is.na(barplot_data$is_evo_force_call_significant),yes = "Non-significant call",no = barplot_data$is_evo_force_call_significant)
barplot_data$type_evo_force <- ifelse(test = is.na(barplot_data$type_evo_force),yes = "Unclear signature",no = barplot_data$type_evo_force)

ggplot(barplot_data, aes(x = type_evo_force, y = n, fill = is_evo_force_call_significant)) +
  geom_col(position = "dodge") +
  facet_grid(is_power_enough_for_pathway_linreg ~ is_power_enough_for_sp_pairs_pathway_based_dist_comparison) +
  xlab("Type of evolutionary force acting at\nthe transcriptomic level")+
  ylab("Number of KEGG pathways (out of 411)")+
  theme(axis.text.x= element_text(angle = 90, hjust = 1,vjust = 0.5), legend.position = "bottom")+
  labs(title = "", fill = "Significance of the evolutionary call")
ggsave(filename = "Barplot_evo_force_type_vs_power_bins_KEGG_pathways.svg", path=output_workspace, width = 17.4/2.54, height = 17.4/2.54,device = svg)
ggsave(filename = "Barplot_evo_force_type_vs_power_bins_KEGG_pathways.png", path=output_workspace, width = 17.4/2.54, height = 17.4/2.54,device = png,dpi = 300)

#Illustrate the frequency of evolutionary forces across power bins using a ggplot pie chart
  #Make grid of piecharts
plot_data <- df_ko_pathways_evo_force %>%
  count(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison, most_likely_and_signif_evo_force) %>%
  group_by(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison) %>%
  mutate(percent = n / sum(n))

ggplot(plot_data, aes(x = "", y = percent, fill = most_likely_and_signif_evo_force)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  facet_grid(is_power_enough_for_pathway_linreg ~ is_power_enough_for_sp_pairs_pathway_based_dist_comparison) +
  labs(title = "", fill = "Evolutionary call\nfor KEGG pathways")
ggsave(filename = "Piechart_evo_calls_vs_power_bins.svg", path=output_workspace, width = 17.4/2.54, height = 13.2/2.54,device = svg)

#make it as a barplot
barplot_data <- df_ko_pathways_evo_force %>%
  count(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison, is_evo_force_call_significant, most_likely_evo_force) %>%
  group_by(is_power_enough_for_pathway_linreg, is_power_enough_for_sp_pairs_pathway_based_dist_comparison,is_evo_force_call_significant) %>%
  mutate(percent = n / sum(n))

barplot_data$is_evo_force_call_significant <- ifelse(test = is.na(barplot_data$is_evo_force_call_significant),yes = "Non-significant call",no = barplot_data$is_evo_force_call_significant)
barplot_data$most_likely_evo_force <- ifelse(test = is.na(barplot_data$most_likely_evo_force),yes = "Unclear signature",no = barplot_data$most_likely_evo_force)

ggplot(barplot_data, aes(x = most_likely_evo_force, y = n, fill = is_evo_force_call_significant)) +
  geom_col(position = "dodge") +
  facet_grid(is_power_enough_for_pathway_linreg ~ is_power_enough_for_sp_pairs_pathway_based_dist_comparison) +
  xlab("Evolutionary call")+
  ylab("Number of KEGG pathways (out of 411)")+
  theme(axis.text.x= element_text(angle = 90, hjust = 1,vjust = 0.5), legend.position = "bottom")+
  labs(title = "", fill = "Significance of the evolutionary call")
ggsave(filename = "Barplot_evo_force_vs_power_bins_KEGG_pathways.svg", path=output_workspace, width = 17.4/2.54, height = 17.4/2.54,device = svg)
ggsave(filename = "Barplot_evo_force_vs_power_bins_KEGG_pathways.png", path=output_workspace, width = 17.4/2.54, height = 17.4/2.54,device = png,dpi = 300)

#Summarize evolutionary calls
df_ko_pathways_evo_force$most_likely_evo_force <- ifelse(test = is.na(df_ko_pathways_evo_force$most_likely_evo_force),yes = "Unclear signature",no = df_ko_pathways_evo_force$most_likely_evo_force)
df_ko_pathways_evo_force$is_evo_force_call_supported_by_correlation_analysis <- ifelse(test=df_ko_pathways_evo_force$most_likely_evo_force == "Unclear signature", yes= F, no=T)
df_ko_pathways_evo_force$is_evo_force_call_supported_by_distance_analysis <- ifelse(test = df_ko_pathways_evo_force$most_likely_evo_force == "Unclear signature",yes = ifelse(test = (df_ko_pathways_evo_force$pval_evo_force_call >= 0.05)||is.na(df_ko_pathways_evo_force$pval_evo_force_call),yes = T,no = F),no = ifelse(test = df_ko_pathways_evo_force$pval_evo_force_call < 0.05,yes = T,no = F))
df_summary_ko_pathways_evo_calls <-df_ko_pathways_evo_force %>%
  count(most_likely_evo_force,type_evo_force,is_evo_force_call_supported_by_correlation_analysis, is_evo_force_call_supported_by_distance_analysis) %>%
  mutate(percent = n / sum(n))

#save full evo call table and summary table
write.table(x=df_ko_pathways_evo_force,file = paste0(output_workspace,"YPD_table_df_ko_pathways_evo_force.tsv"),sep = "\t",na = "NA",row.names = F,col.names = T,quote = F)
write.table(x=df_summary_ko_pathways_evo_calls,file = paste0(output_workspace,"YPD_table_df_summary_ko_pathways_evo_calls.tsv"),sep = "\t",na = "NA",row.names = F,col.names = T,quote = F)

################Investigate the structure of genes and strains based on presence/absence data + expression data
  #16 Yeasts species genes CLUSTERING BASED ON PRESENCE/ABSENCE across strains(are some orthogroups acquired/lost together and community gene content or geographic isolation effect or mobile elements?)
    #Gower distance for rows
row_dist_sp_accessory_orthogroup_pres_abs_dist <- vegdist(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[,v_lst_accessory_Orthogroups_YPD], method = "gower",na.rm=T)
    #Gower distance for columns (transpose the matrix first)
col_dist_sp_accessory_orthogroup_pres_abs_dist <- vegdist(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[,v_lst_accessory_Orthogroups_YPD]), method = "gower",na.rm=T)
    #Ward clustering for rows
row_clust_sp_accessory_orthogroup_pres_abs_dist <- hclust(row_dist_sp_accessory_orthogroup_pres_abs_dist, method = "ward.D2") #Use ward.D2 for Ward's method
    #Ward clustering for columns
col_clust_sp_accessory_orthogroup_pres_abs_dist <- hclust(col_dist_sp_accessory_orthogroup_pres_abs_dist, method = "ward.D2")

#Draw silhouette score profile
dist_matrix_accessory_orthogroups <- col_dist_sp_accessory_orthogroup_pres_abs_dist #Calculate distance matrix
hclust_result_accessory_orthogroups <- col_clust_sp_accessory_orthogroup_pres_abs_dist #Perform hierarchical clustering
k_range_accessory_orthogroups <- seq(10,3000,10)
cluster_assignments_accessory_orthogroups <- list()
for (k in k_range_accessory_orthogroups) {
  cluster_assignments_accessory_orthogroups[[as.character(k)]] <- cutree(hclust_result_accessory_orthogroups, k = k)
}
silhouette_scores_accessory_orthogroups <- numeric(length(k_range_accessory_orthogroups))
for (i in seq_along(k_range_accessory_orthogroups)) {
  k <- k_range_accessory_orthogroups[i]
  sil_obj <- cluster::silhouette(cluster_assignments_accessory_orthogroups[[as.character(k)]], dist_matrix_accessory_orthogroups)
  silhouette_scores_accessory_orthogroups[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_accessory_orthogroups_pres_abs_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_accessory_orthogroups, silhouette_scores_accessory_orthogroups, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
#Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_acc_orthogroups <- k_range_accessory_orthogroups[silhouette_scores_accessory_orthogroups==max(silhouette_scores_accessory_orthogroups)]
clusters_vector_acc_orthogroups <- cutree(hclust_result_accessory_orthogroups, k = nb_clusters_acc_orthogroups)
    #Test the significance of the clusters
      #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_acc_orthogroups, dist_matrix_accessory_orthogroups)[, 3])
      #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_acc_orthogroups)
  mean(silhouette(shuffled_clusters, dist_matrix_accessory_orthogroups)[, 3])
})
      #Calculate p-value
p_value_signif_clusters_acc_orthogroups <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_acc_orthogroups


  ####~1000 Scer gene CLUSTERING BASED ON PRESENCE/ABSENCE across strains(are some genes acquired/lost together?)
    #Scer vectors of accessory and core genes
v_lst_Scer_accessory_genes_with_expr_data <- rownames(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data))[rowSums(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data),na.rm = T)/ncol(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data))!=1]
v_lst_Scer_core_genes_with_expr_data <- rownames(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data))[rowSums(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data),na.rm = T)/ncol(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data))==1]
mtx_log2_1p_Scer_gene_expr_TPM_YPD <- log2(mtx_Scer_expression_for_genes_with_pres_abs_data+1)

      #Scer gene expression abundance and dispersion figure
v_abundance_accessory_genes_Scer_YPD <- apply(X = mtx_log2_1p_Scer_gene_expr_TPM_YPD[,v_lst_Scer_accessory_genes_with_expr_data],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_abundance_core_genes_Scer_YPD <- apply(X = mtx_log2_1p_Scer_gene_expr_TPM_YPD[,v_lst_Scer_core_genes_with_expr_data],MARGIN = 2,FUN = function(x) mean(x,na.rm=T))
v_dispersion_accessory_genes_Scer_YPD <- apply(X = mtx_log2_1p_Scer_gene_expr_TPM_YPD[,v_lst_Scer_accessory_genes_with_expr_data],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))
v_dispersion_core_genes_Scer_YPD <- apply(X = mtx_log2_1p_Scer_gene_expr_TPM_YPD[,v_lst_Scer_core_genes_with_expr_data],MARGIN = 2,FUN = function(x) ie2misc::madstat(x,na.rm=T))

df_expr_abundance_pangenome_Yeast_spcs_vs_Scer <- data.frame(Gene=c(names(v_abundance_accessory_genes_Scer_YPD),names(v_dispersion_accessory_genes_Scer_YPD),names(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),names(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),names(v_abundance_core_genes_Scer_YPD),names(v_dispersion_core_genes_Scer_YPD),names(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),names(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD)),
                                                        Dataset=rep(c("Baker's yeast strains","16 yeast species","Baker's yeast strains","16 yeast species"),times=c(length(v_abundance_accessory_genes_Scer_YPD)+length(v_dispersion_accessory_genes_Scer_YPD),length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_core_genes_Scer_YPD)+length(v_dispersion_core_genes_Scer_YPD),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD) )),
                                                        Metric=rep(c("Abundance (Expression level)","Dispersion (Mean Absolute Deviation)","Abundance (Expression level)","Dispersion (Mean Absolute Deviation)","Abundance (Expression level)","Dispersion (Mean Absolute Deviation)","Abundance (Expression level)","Dispersion (Mean Absolute Deviation)"),times=c(length(v_abundance_accessory_genes_Scer_YPD),length(v_dispersion_accessory_genes_Scer_YPD),length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_core_genes_Scer_YPD),length(v_dispersion_core_genes_Scer_YPD),length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD) )),
                                                        Pangenomic_component=rep(c("Accessory","Core"),times=c(length(v_abundance_accessory_genes_Scer_YPD)+length(v_dispersion_accessory_genes_Scer_YPD)+length(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),length(v_abundance_core_genes_Scer_YPD)+length(v_dispersion_core_genes_Scer_YPD)+length(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD)+length(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD) )),
                                                        Value=c(unname(v_abundance_accessory_genes_Scer_YPD),unname(v_dispersion_accessory_genes_Scer_YPD),unname(v_abundance_accessory_Orthogroups_Yeast_21_spcs_YPD),unname(v_dispersion_accessory_Orthogroups_Yeast_21_spcs_YPD),unname(v_abundance_core_genes_Scer_YPD),unname(v_dispersion_core_genes_Scer_YPD),unname(v_abundance_core_Orthogroups_Yeast_21_spcs_YPD),unname(v_dispersion_core_Orthogroups_Yeast_21_spcs_YPD) ),stringsAsFactors = F)
df_expr_abundance_pangenome_Yeast_spcs_vs_Scer$Dataset <- factor(df_expr_abundance_pangenome_Yeast_spcs_vs_Scer$Dataset,levels=c("Baker's yeast strains","16 yeast species"))
df_expr_abundance_pangenome_Yeast_spcs_vs_Scer$Pangenomic_component <- factor(df_expr_abundance_pangenome_Yeast_spcs_vs_Scer$Pangenomic_component,levels=c("Core","Accessory"))
ggplot(data = df_expr_abundance_pangenome_Yeast_spcs_vs_Scer,aes(x=as.factor(Pangenomic_component),y = Value)) + geom_violin(fill="grey") + geom_boxplot(width=0.075) + xlab("Pangenomic component") + ylab("Value") + theme_bw() + theme(axis.title = element_text(size=12),axis.text = element_text(size=12),strip.text = element_text(size = 14)) + stat_compare_means(method = "wilcox") + facet_grid(Metric~Dataset,scales="free_y")
ggsave(filename = "Expression_Abundance_AND_Dispersion_Accessory_vs_Core_Orthogroup_in_Scer_vs_16YeastSpecies.png", path=output_workspace, width = 25, height = 20, units = "cm",dpi = 1200)
ggsave(filename = "Expression_Abundance_AND_Dispersion_Accessory_vs_Core_Orthogroup_in_Scer_vs_16YeastSpecies.svg", path=output_workspace, width = 25/2.54, height = 20/2.54, device=svg)

    #Gower distance for rows
row_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist <- vegdist(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data[,v_lst_Scer_accessory_genes_with_expr_data], method = "gower",na.rm=T)
    #Gower distance for columns (transpose the matrix first)
      #extract genes with non-problematic presence/absence profile (can be compared to other genes because NA or 0 are less frequent than average; This filter preserves most of the genes as the problematic ones have outliers number of NAs and 0s)
col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist <- vegdist(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data[,v_lst_Scer_accessory_genes_with_expr_data]), method = "gower",na.rm=T)
v_extract_non_problematic_accessory_Genes <- rowSums(is.na(as.matrix(col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist)))<mean(rowSums(is.na(as.matrix(col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist))))
col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist <- as.matrix(col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist)[v_extract_non_problematic_accessory_Genes,v_extract_non_problematic_accessory_Genes]
    #Ward clustering for rows
row_clust_Scer_accessory_genes_with_expr_data_pres_abs_dist <- hclust(row_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist, method = "ward.D2") #Use ward.D2 for Ward's method
    #Ward clustering for columns
col_clust_Scer_accessory_genes_with_expr_data_pres_abs_dist <- hclust(as.dist(col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist), method = "ward.D2")

    #Draw silhouette score profile
dist_matrix_Scer_accessory_genes <- col_dist_Scer_accessory_genes_with_expr_data_pres_abs_dist #Calculate distance matrix
hclust_result_Scer_accessory_genes <- col_clust_Scer_accessory_genes_with_expr_data_pres_abs_dist #Perform hierarchical clustering
k_range_Scer_accessory_genes <- seq(5,1200,5)
cluster_assignments_Scer_accessory_genes <- list()
for (k in k_range_Scer_accessory_genes) {
  cluster_assignments_Scer_accessory_genes[[as.character(k)]] <- cutree(hclust_result_Scer_accessory_genes, k = k)
}
silhouette_scores_Scer_accessory_genes_pres_abs <- numeric(length(k_range_Scer_accessory_genes))
for (i in seq_along(k_range_Scer_accessory_genes)) {
  k <- k_range_Scer_accessory_genes[i]
  sil_obj <- cluster::silhouette(cluster_assignments_Scer_accessory_genes[[as.character(k)]], dist_matrix_Scer_accessory_genes)
  silhouette_scores_Scer_accessory_genes_pres_abs[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_Scer_accessory_genes_pres_abs_pres_abs_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_Scer_accessory_genes, silhouette_scores_Scer_accessory_genes_pres_abs, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()

    #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_Scer_acc_genes_pres_abs <- k_range_Scer_accessory_genes[silhouette_scores_Scer_accessory_genes_pres_abs==max(silhouette_scores_Scer_accessory_genes_pres_abs)]
clusters_vector_Scer_acc_genes_pres_abs <- cutree(hclust_result_Scer_accessory_genes, k = nb_clusters_Scer_acc_genes_pres_abs)
    #Test the significance of the clusters
    #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_Scer_acc_genes_pres_abs, dist_matrix_Scer_accessory_genes)[, 3])
    #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_Scer_acc_genes_pres_abs)
  mean(silhouette(shuffled_clusters, dist_matrix_Scer_accessory_genes)[, 3])
})
    #Calculate p-value
p_value_signif_clusters_Scer_acc_genes_based_on_pres_abs_only <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_Scer_acc_genes_based_on_pres_abs_only

  #In S. cerevisiae, assess the effect of population structure on expression  (Strain expression cluster vs Strain location of origin)
    #Import ~1000 Scer strains metadata (including location of origin)
df_Scer_strains_metadata <- read.csv2(file = paste0(output_workspace,"Scer_strains_metadata.tsv"),sep = "\t",header = T,quote="",stringsAsFactors = FALSE)
      #keep country only when possible
for (i in 1:nrow(df_Scer_strains_metadata)){
  str_orig_loc <- df_Scer_strains_metadata$Geographical.origins[i]
  v_sep_loc <- strsplit(x = str_orig_loc,split = ", ",fixed = T)[[1]]
  df_Scer_strains_metadata$Geographical.origins[i] <- gsub(pattern = ",",replacement = "",x = gsub(pattern = " ",replacement = "",x = v_sep_loc[length(v_sep_loc)],fixed = T),fixed = T)
}

v_Scer_strain_to_geographical_origin <- df_Scer_strains_metadata$Geographical.origins
names(v_Scer_strain_to_geographical_origin) <- df_Scer_strains_metadata$Standardized.name

    #Visualize ~1000 Scer strains expression PCA of the genes, mapping back species + color species by continent or country
df_Scer_regular_gene_expr_pca$geographical_origin <- v_Scer_strain_to_geographical_origin[df_Scer_regular_gene_expr_pca$species]
df_Scer_regular_gene_expr_pca %>% 
  ggplot(aes(x = PC1, y = PC2, col=geographical_origin)) +
  labs(x = "PC1", y = "PC2")+
  geom_point(size = 2, stroke = 2) + 
  theme(legend.position = "none") +
  ggrepel::geom_text_repel(aes(label = mixture)) #geom_text
ggsave(filename = "PCA_Scer_regular_gene_expr.png", path=output_workspace, width = 35/2.54, height = 25/2.54)

    ####create clusters of accessory genes based on expression data
      #Bray-Curtis distance for rows
row_dist_Scer_accessory_genes_with_expr_data_dist <- vegdist(mtx_Scer_expression_for_genes_with_pres_abs_data[,v_lst_Scer_accessory_genes_with_expr_data], method = "bray",na.rm=T)
      ##Bray-Curtis distance for columns (transpose the matrix first)
col_dist_Scer_accessory_genes_with_expr_data_dist <- vegdist(t(mtx_Scer_expression_for_genes_with_pres_abs_data[,v_lst_Scer_accessory_genes_with_expr_data]), method = "bray",na.rm=T)
          #remove problematic accessory genes (cant be compared with other genes because too rare)
v_extract_non_problematic_accessory_Genes <- rowSums(is.na(as.matrix(col_dist_Scer_accessory_genes_with_expr_data_dist)))==0
col_dist_Scer_accessory_genes_with_expr_data_dist <- as.dist(as.matrix(col_dist_Scer_accessory_genes_with_expr_data_dist)[v_extract_non_problematic_accessory_Genes,v_extract_non_problematic_accessory_Genes])
          #dim(as.matrix(col_dist_Scer_accessory_genes_with_expr_data_dist))
      #Ward clustering for rows
row_clust_Scer_accessory_genes_with_expr_data_dist <- hclust(row_dist_Scer_accessory_genes_with_expr_data_dist, method = "ward.D2") #Use ward.D2 for Ward's method
      #Ward clustering for columns
col_clust_Scer_accessory_genes_with_expr_data_dist <- hclust(col_dist_Scer_accessory_genes_with_expr_data_dist, method = "ward.D2")
      #Draw silhouette score profile
dist_matrix_Scer_accessory_genes <- col_dist_Scer_accessory_genes_with_expr_data_dist #Calculate distance matrix
hclust_result_Scer_accessory_genes <- col_clust_Scer_accessory_genes_with_expr_data_dist #Perform hierarchical clustering
k_range_Scer_accessory_genes <- seq(10,1200,10)
cluster_assignments_Scer_accessory_genes <- list()
for (k in k_range_Scer_accessory_genes) {
  cluster_assignments_Scer_accessory_genes[[as.character(k)]] <- cutree(hclust_result_Scer_accessory_genes, k = k)
}
silhouette_scores_Scer_accessory_genes_expr_data <- numeric(length(k_range_Scer_accessory_genes))
for (i in seq_along(k_range_Scer_accessory_genes)) {
  k <- k_range_Scer_accessory_genes[i]
  sil_obj <- cluster::silhouette(cluster_assignments_Scer_accessory_genes[[as.character(k)]], dist_matrix_Scer_accessory_genes)
  silhouette_scores_Scer_accessory_genes_expr_data[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_Scer_accessory_genes_expr_data_expr_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_Scer_accessory_genes, silhouette_scores_Scer_accessory_genes_expr_data, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
      #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_Scer_acc_genes_expr_data <- k_range_Scer_accessory_genes[silhouette_scores_Scer_accessory_genes_expr_data==max(silhouette_scores_Scer_accessory_genes_expr_data)]
clusters_vector_Scer_acc_genes_expr_data <- cutree(hclust_result_Scer_accessory_genes, k = nb_clusters_Scer_acc_genes_expr_data)
      #Test the significance of the clusters
        #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_Scer_acc_genes_expr_data, dist_matrix_Scer_accessory_genes)[, 3])
      #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_Scer_acc_genes_expr_data)
  mean(silhouette(shuffled_clusters, dist_matrix_Scer_accessory_genes)[, 3])
})
      #Calculate p-value
p_value_signif_clusters_Scer_acc_genes_based_on_expr_data <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_Scer_acc_genes_based_on_expr_data

  ##create clusters of strains based on expression data
        #Draw silhouette score profile
dist_matrix_Scer_strains_based_on_expr_data <- row_dist_Scer_accessory_genes_with_expr_data_dist #Calculate distance matrix
hclust_result_Scer_strains_based_on_expr_data <- row_clust_Scer_accessory_genes_with_expr_data_dist #Perform hierarchical clustering
k_range_Scer_strains <- seq(5,700,5)
cluster_assignments_Scer_strains <- list()
for (k in k_range_Scer_strains) {
  cluster_assignments_Scer_strains[[as.character(k)]] <- cutree(hclust_result_Scer_strains_based_on_expr_data, k = k)
}
silhouette_scores_Scer_strains <- numeric(length(k_range_Scer_strains))
for (i in seq_along(k_range_Scer_strains)) {
  k <- k_range_Scer_strains[i]
  sil_obj <- cluster::silhouette(cluster_assignments_Scer_strains[[as.character(k)]], dist_matrix_Scer_strains_based_on_expr_data)
  silhouette_scores_Scer_strains[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_Scer_strains_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_Scer_strains, silhouette_scores_Scer_strains, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
        #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_Scer_strains_based_on_expr_data <- k_range_Scer_strains[silhouette_scores_Scer_strains==max(silhouette_scores_Scer_strains)]
clusters_vector_Scer_strains_based_on_expr_data <- cutree(hclust_result_Scer_strains_based_on_expr_data, k = nb_clusters_Scer_strains_based_on_expr_data)
        #Test the significance of the clusters
          #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_Scer_strains_based_on_expr_data, dist_matrix_Scer_strains_based_on_expr_data)[, 3])
        #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_Scer_strains_based_on_expr_data)
  mean(silhouette(shuffled_clusters, dist_matrix_Scer_strains_based_on_expr_data)[, 3])
})
        #Calculate p-value
p_value_signif_clusters_Scer_strains_based_on_expr_data <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_Scer_strains_based_on_expr_data
    #In S. cerevisiae, test the association between strain cluster and strain geographical origin
v_strains_location <- v_Scer_strain_to_geographical_origin[rownames(as.matrix(row_dist_Scer_accessory_genes_with_expr_data_dist))]
pval_association_Scer_strain_Expr_cluster_vs_geographical_origin <- chisq.test(x=clusters_vector_Scer_strains_based_on_expr_data[!is.na(v_strains_location)],y = v_strains_location[!is.na(v_strains_location)])
pval_association_Scer_strain_Expr_cluster_vs_geographical_origin
lsr::cramersV(table(clusters_vector_Scer_strains_based_on_expr_data[!is.na(v_strains_location)],v_strains_location[!is.na(v_strains_location)]))
    #strain clusters vs strain geographical origin
table(clusters_vector_Scer_strains_based_on_expr_data[!is.na(v_strains_location)], v_strains_location[!is.na(v_strains_location)])
#View(subset(as.data.frame(table(clusters_vector_Scer_strains_based_on_expr_data[!is.na(v_strains_location)], v_strains_location[!is.na(v_strains_location)])), Freq>0))
##create clusters of genes (WHOLE PANGENOME NOT JUST ACCESSORY) based on expression data
        #Bray-Curtis distance for rows
row_dist_Scer_ALL_genes_with_expr_data_dist <- vegdist(mtx_Scer_expression_for_genes_with_pres_abs_data, method = "bray",na.rm=T)
        ##Bray-Curtis distance for columns (transpose the matrix first)
col_dist_Scer_ALL_genes_with_expr_data_dist <- vegdist(t(mtx_Scer_expression_for_genes_with_pres_abs_data), method = "bray",na.rm=T)
        #remove problematic genes (cant be compared with other genes because too rare)
v_extract_non_problematic_Genes <- rowSums(is.na(as.matrix(col_dist_Scer_ALL_genes_with_expr_data_dist)))==0
col_dist_Scer_ALL_genes_with_expr_data_dist <- as.dist(as.matrix(col_dist_Scer_ALL_genes_with_expr_data_dist)[v_extract_non_problematic_Genes,v_extract_non_problematic_Genes])
#dim(as.matrix(col_dist_Scer_ALL_genes_with_expr_data_dist))
        #Ward clustering for rows
row_clust_Scer_ALL_genes_with_expr_data_dist <- hclust(row_dist_Scer_ALL_genes_with_expr_data_dist, method = "ward.D2") #Use ward.D2 for Ward's method
        #Ward clustering for columns
col_clust_Scer_ALL_genes_with_expr_data_dist <- hclust(col_dist_Scer_ALL_genes_with_expr_data_dist, method = "ward.D2")
        #Draw silhouette score profile
dist_matrix_Scer_ALL_genes <- col_dist_Scer_ALL_genes_with_expr_data_dist #Calculate distance matrix
hclust_result_Scer_ALL_genes <- col_clust_Scer_ALL_genes_with_expr_data_dist #Perform hierarchical clustering
k_range_Scer_ALL_genes <- seq(2,1200,3)
cluster_assignments_Scer_ALL_genes <- list()
for (k in k_range_Scer_ALL_genes) {
  cluster_assignments_Scer_ALL_genes[[as.character(k)]] <- cutree(hclust_result_Scer_ALL_genes, k = k)
}
silhouette_scores_Scer_ALL_genes_expr_data <- numeric(length(k_range_Scer_ALL_genes))
for (i in seq_along(k_range_Scer_ALL_genes)) {
  k <- k_range_Scer_ALL_genes[i]
  sil_obj <- cluster::silhouette(cluster_assignments_Scer_ALL_genes[[as.character(k)]], dist_matrix_Scer_ALL_genes)
  silhouette_scores_Scer_ALL_genes_expr_data[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_Scer_ALL_genes_expr_data_expr_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_Scer_ALL_genes, silhouette_scores_Scer_ALL_genes_expr_data, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
        #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_Scer_ALL_genes_expr_data <- k_range_Scer_ALL_genes[silhouette_scores_Scer_ALL_genes_expr_data==max(silhouette_scores_Scer_ALL_genes_expr_data)]
clusters_vector_Scer_ALL_genes_expr_data <- cutree(hclust_result_Scer_ALL_genes, k = nb_clusters_Scer_ALL_genes_expr_data)
        #Test the significance of the clusters
          #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_Scer_ALL_genes_expr_data, dist_matrix_Scer_ALL_genes)[, 3])
        #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_Scer_ALL_genes_expr_data)
  mean(silhouette(shuffled_clusters, dist_matrix_Scer_ALL_genes)[, 3])
})
        #Calculate p-value
p_value_signif_clusters_Scer_ALL_genes_based_on_expr_data <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_Scer_ALL_genes_based_on_expr_data
        #In S. cerevisiae, test the association between expression clusters of genes and gene class (accessory/core genes)
v_is_core_tested_genes <- colnames(as.matrix(col_dist_Scer_ALL_genes_with_expr_data_dist))%in%v_lst_Scer_core_genes_with_expr_data
pval_association_Scer_ALL_genes_Expr_cluster_vs_gene_class <- chisq.test(x = clusters_vector_Scer_ALL_genes_expr_data,y = v_is_core_tested_genes)$p.value
lsr::cramersV(table(clusters_vector_Scer_ALL_genes_expr_data,v_is_core_tested_genes))
          #expression clusters vs is_the_gene_core
table(clusters_vector_Scer_ALL_genes_expr_data, v_is_core_tested_genes)
          #Heatmap expression clusters vs is_the_gene_core
row_dist <- row_dist_Scer_ALL_genes_with_expr_data_dist
col_dist <- col_dist_Scer_ALL_genes_with_expr_data_dist
row_clust <- hclust(row_dist, method = "ward.D2") #Use ward.D2 for Ward's method
col_clust <- hclust(col_dist, method = "ward.D2")
    #assign colors to strains locations of origin and pangenome components
v_strains_location_cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(sort(unique(v_strains_location))))
names(v_strains_location_cols) <- sort(unique(v_strains_location))
v_cols_rows_strains <- unname(v_strains_location_cols[v_strains_location])
v_cols_rows_strains[is.na(v_cols_rows_strains)] <- "white"
png(filename = paste0(output_workspace,"Clustered_Heatmap_ALL_Genes_expr_and_pangenome.png"),width = 25,height = 25,units = "cm",res=300)
heatmap.2(log2(mtx_Scer_expression_for_genes_with_pres_abs_data+1)[,v_extract_non_problematic_Genes],
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust),
          Colv = as.dendrogram(col_clust),
          dendrogram = "both",
          trace = "none",
          main = NULL,
          xlab = "Genes",
          ylab="Baker's yeast strains",
          labRow=NA,
          labCol=NA,
          key.par = list(cex = 0.45),
          key.title = NA,
          key.xlab = "log2(TPM+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          margins = c(2, 2),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2),
          RowSideColors=v_cols_rows_strains,
          ColSideColors=ifelse(test = v_is_core_tested_genes,yes = "red3",no="skyblue")
          )
dev.off()
svg(filename = paste0(output_workspace,"Clustered_Heatmap_ALL_Genes_expr_and_pangenome.svg"),width = 25,height = 25,units = "cm")
heatmap.2(log2(mtx_Scer_expression_for_genes_with_pres_abs_data+1)[,v_extract_non_problematic_Genes],
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust),
          Colv = as.dendrogram(col_clust),
          dendrogram = "both",
          trace = "none",
          main = NULL,
          xlab = "Genes",
          ylab="Baker's yeast strains",
          labRow=NA,
          labCol=NA,
          key.par = list(cex = 0.45),
          key.title = NA,
          key.xlab = "log2(TPM+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          margins = c(2, 2),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2),
          RowSideColors=v_cols_rows_strains,
          ColSideColors=ifelse(test = v_is_core_tested_genes,yes = "red3",no="skyblue")
)
dev.off()

      #S. cerevisiae strains' presence absence heatmap
        #Gower distance for rows
row_dist_Scer_all_genes_with_expr_data_pres_abs_dist <- vegdist(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data, method = "gower",na.rm=T)
        #Gower distance for columns (transpose the matrix first)
          #extract genes with non-problematic presence/absence profile (can be compared to other genes because NA or 0 are less frequent than average; This filter preserves most of the genes as the problematic ones have outliers number of NAs and 0s)
col_dist_Scer_all_genes_with_expr_data_pres_abs_dist <- vegdist(t(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data), method = "gower",na.rm=T)
v_extract_non_problematic_Genes <- rowSums(is.na(as.matrix(col_dist_Scer_all_genes_with_expr_data_pres_abs_dist)))<mean(rowSums(is.na(as.matrix(col_dist_Scer_all_genes_with_expr_data_pres_abs_dist))))
col_dist_Scer_all_genes_with_expr_data_pres_abs_dist <- as.matrix(col_dist_Scer_all_genes_with_expr_data_pres_abs_dist)[v_extract_non_problematic_Genes,v_extract_non_problematic_Genes]
          #Ward clustering for rows
row_clust_Scer_all_genes_with_expr_data_pres_abs_dist <- hclust(row_dist_Scer_all_genes_with_expr_data_pres_abs_dist, method = "ward.D2") #Use ward.D2 for Ward's method
          #Ward clustering for columns
col_clust_Scer_all_genes_with_expr_data_pres_abs_dist <- hclust(as.dist(col_dist_Scer_all_genes_with_expr_data_pres_abs_dist), method = "ward.D2")
png(filename = paste0(output_workspace,"Clustered_Heatmap_Pres_Abs_ALL_Genes.png"),width = 25,height = 25,units = "cm",res=300)
heatmap.2(as.matrix(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data)[,v_extract_non_problematic_Genes],
          distfun = function(x) vegdist(x, method = "gower"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust_Scer_all_genes_with_expr_data_pres_abs_dist),
          Colv = as.dendrogram(col_clust_Scer_all_genes_with_expr_data_pres_abs_dist),
          dendrogram = "both",
          trace = "none",
          main = NULL,
          xlab = "Genes",
          ylab="Baker's yeast strains",
          labRow=NA,
          labCol=NA,
          key.par = list(cex = 0.45),
          key.title = NA,
          key.xlab = "Presence",
          key.ylab = "Count",
          col= c("white","black"),
          na.color = "gray", 
          margins = c(2, 2),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2)
)
dev.off()

    #List Scer genes highly expressed in specific strains
the_mtx <- as.matrix(mtx_Scer_pres_abs_filtered_for_genes_with_expr_data)[,v_extract_non_problematic_Genes]
global_mean <- mean(the_mtx, na.rm = TRUE)
n_species <- nrow(the_mtx)
pct_above_mean <- colMeans(the_mtx > global_mean, na.rm = TRUE)
pct_is_zero <- colMeans(the_mtx == 0, na.rm = TRUE)
Scer_genes_highly_expressed_in_specific_strains <- colnames(the_mtx)[
  pct_above_mean >= 0.10 & 
    pct_above_mean <= 0.50 & 
    pct_is_zero >= 0.50
]
Scer_genes_highly_expressed_in_specific_strains <- Scer_genes_highly_expressed_in_specific_strains[!is.na(Scer_genes_highly_expressed_in_specific_strains)]
write.table(x=Scer_genes_highly_expressed_in_specific_strains,file = paste0(output_workspace,"List_Scer_genes_highly_expressed_in_specific_strains.tsv"),sep = "\t",na = "NA",row.names = F,col.names = F,quote = F)

    ###Now analyse clusters for the 16 YPD yeast species with the orthogroups
      #Clustered heatmap for orthogroup expression vs pangenome structure colsidecolorbar
      ##create clusters of orthogroups (WHOLE PANGENOME NOT JUST ACCESSORY) based on expression data
        #Bray-Curtis distance for rows
row_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist <- vegdist(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD, method = "bray",na.rm=T)
        ##Bray-Curtis distance for columns (transpose the matrix first)
col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist <- vegdist(t(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD), method = "bray",na.rm=T)
        #remove problematic orthogroups (cant be compared with other orthogroups because too rare)
v_extract_non_problematic_orthogroups <- rowSums(is.na(as.matrix(col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist)))==0
col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist <- as.dist(as.matrix(col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist)[v_extract_non_problematic_orthogroups,v_extract_non_problematic_orthogroups])
        #dim(as.matrix(col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist))
        #Ward clustering for rows
row_clust_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist <- hclust(row_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist, method = "ward.D2") #Use ward.D2 for Ward's method
        #Ward clustering for columns
col_clust_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist <- hclust(col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist, method = "ward.D2")
        #Draw silhouette score profile
dist_matrix_16_YPD_spcs_ALL_orthogroups <- col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist #Calculate distance matrix
hclust_result_16_YPD_spcs_ALL_orthogroups <- col_clust_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist #Perform hierarchical clustering
k_range_16_YPD_spcs_ALL_orthogroups <- seq(2,1200,3)
cluster_assignments_16_YPD_spcs_ALL_orthogroups <- list()
for (k in k_range_16_YPD_spcs_ALL_orthogroups) {
  cluster_assignments_16_YPD_spcs_ALL_orthogroups[[as.character(k)]] <- cutree(hclust_result_16_YPD_spcs_ALL_orthogroups, k = k)
}
silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data <- numeric(length(k_range_16_YPD_spcs_ALL_orthogroups))
for (i in seq_along(k_range_16_YPD_spcs_ALL_orthogroups)) {
  k <- k_range_16_YPD_spcs_ALL_orthogroups[i]
  sil_obj <- cluster::silhouette(cluster_assignments_16_YPD_spcs_ALL_orthogroups[[as.character(k)]], dist_matrix_16_YPD_spcs_ALL_orthogroups)
  silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data[i] <- mean(sil_obj[, "sil_width"])
}
png(filename = paste0(output_workspace,"silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data_expr_clustering.png"),width = 25,height = 25,units = "cm",res = 300)
plot(k_range_16_YPD_spcs_ALL_orthogroups, silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data, type = "o",
     xlab = "Number of Clusters (k)",
     ylab = "Average Silhouette Width",
     main = "Average Silhouette Width vs. Number of Clusters",
     pch = 19, col = "blue")
dev.off()
      #Extract row clusters into k groups that optimizes the silhouette score
nb_clusters_16_YPD_spcs_ALL_orthogroups_expr_data <- k_range_16_YPD_spcs_ALL_orthogroups[silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data==max(silhouette_scores_16_YPD_spcs_ALL_orthogroups_expr_data)]
clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data <- cutree(hclust_result_16_YPD_spcs_ALL_orthogroups, k = nb_clusters_16_YPD_spcs_ALL_orthogroups_expr_data)
      #Test the significance of the clusters
        #Calculate the observed silhouette score
obs_sil <- mean(silhouette(clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data, dist_matrix_16_YPD_spcs_ALL_orthogroups)[, 3])
        #Get the permuted clusters' silhouette scores
perm_sils <- replicate(nb_permutations, {
  shuffled_clusters <- sample(clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data)
  mean(silhouette(shuffled_clusters, dist_matrix_16_YPD_spcs_ALL_orthogroups)[, 3])
})
        #Calculate p-value
p_value_signif_clusters_16_YPD_spcs_ALL_orthogroups_based_on_expr_data <- sum(perm_sils >= obs_sil) / nb_permutations
p_value_signif_clusters_16_YPD_spcs_ALL_orthogroups_based_on_expr_data
      #In S. cerevisiae, test the association between expression clusters of orthogroups and orthogroup class (accessory/core orthogroups)
v_is_core_tested_orthogroups <- colnames(as.matrix(col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist))%in%v_lst_core_Orthogroups_YPD
pval_association_16_YPD_spcs_ALL_orthogroups_Expr_cluster_vs_orthogroup_class <- chisq.test(x = clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data,y = v_is_core_tested_orthogroups)$p.value
lsr::cramersV(table(clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data,v_is_core_tested_orthogroups))
      #expression clusters vs is_the_orthogroup_core
table(clusters_vector_16_YPD_spcs_ALL_orthogroups_expr_data, v_is_core_tested_orthogroups)
      #Heatmap expression clusters vs is_the_orthogroup_core
row_dist <- row_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist
col_dist <- col_dist_16_YPD_spcs_ALL_orthogroups_with_expr_data_dist
row_clust <- hclust(row_dist, method = "ward.D2") #Use ward.D2 for Ward's method
col_clust <- hclust(col_dist, method = "ward.D2")
      #clustered heatmap with yeast subphylum pangenome
svg(filename = paste0(output_workspace,"Clustered_Heatmap_ALL_Orthorgroups_expr_and_pangenome.svg"),width = 25/2.54,height = 25/2.54)
heatmap.2(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups],
          distfun = function(x) vegdist(x, method = "bray"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust),
          Colv = as.dendrogram(col_clust),
          dendrogram = "both",
          trace = "none",
          main = NULL,
          xlab = "Orthogroup",
          ylab="Yeast species",
          labCol=NA,
          cexRow=1,
          key.par = list(cex = 0.45),
          key.title = NA,
          key.xlab = "log2(FPKM+1)",
          key.ylab = "Count",
          col= brewer.pal(n = 9, name = "YlOrRd"),
          margins = c(2, 15),
          lwid = c(1.2, 8),
          lhei = c(0.5, 2),
          RowSideColors = v_clade_color[v_order_almost_21_spcs_YPD[rownames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups])]],
          ColSideColors=ifelse(test = v_is_core_tested_orthogroups,yes = "red3",no="skyblue")
)
dev.off()
        #16 YPD species' orthogroups presence absence heatmap
          #Gower distance for rows
row_dist_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist <- vegdist(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),], method = "gower",na.rm=T)
          #Gower distance for columns (transpose the matrix first)
            #extract orthogroups with non-problematic presence/absence profile (can be compared to other orthogroups because NA or 0 are less frequent than average; This filter preserves most of the orthogroups as the problematic ones have outliers number of NAs and 0s)
col_dist_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist <- vegdist(t(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),]), method = "gower",na.rm=T)
          #Ward clustering for rows
row_clust_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist <- hclust(row_dist_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist, method = "ward.D2") #Use ward.D2 for Ward's method
          #Ward clustering for columns
col_clust_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist <- hclust(as.dist(col_dist_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist), method = "ward.D2")
svg(filename = paste0(output_workspace,"Clustered_Heatmap_Pres_Abs_ALL_orthogroups.svg"),width = 25/2.54,height = 25/2.54)
heatmap.2(as.matrix(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),]),
          distfun = function(x) vegdist(x, method = "gower"), #Specify Bray-Curtis for heatmap.2
          hclustfun = function(x) hclust(x, method = "ward.D2"), #Specify Ward's method
          Rowv = as.dendrogram(row_clust_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist),
          Colv = as.dendrogram(col_clust_16_YPD_spcs_all_orthogroups_with_expr_data_pres_abs_dist),
          dendrogram = "both",
          trace = "none",
          main = NULL,
          xlab = "Orthogroup",
          ylab="Yeast species",
          labCol=NA,
          cexRow=1,
          key.par = list(cex = 0.45),
          key.title = NA,
          key.xlab = "Presence",
          key.ylab = "Count",
          col= c("white","black"),
          na.color = "gray", 
          margins = c(2, 15),
          lmat = rbind(c(0, 5, 4, 0, 0), c(0, 3, 2, 1, 0)), 
          lhei = c(1.5, 4), 
          lwid = c(1, 1, 5, 0.2, 1),
          RowSideColors = v_clade_color[v_order_almost_21_spcs_YPD[rownames(as.matrix(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),]))]],
)
dev.off()

    #Piechart pangenome structure in Scer vs 16 yeast species
      #S. cerevisiae
proportions_Scer_pangenome <- c(proportion_accessory_genes_in_Scer_pangenome,proportion_core_genes_in_Scer_pangenome)
names(proportions_Scer_pangenome) <- c("Accessory\nGenome", "Core\nGenome")
#svg(filename = paste0(output_workspace,"Scer_pangenome_structure.svg"),width = 7.5/2.54,height = 7.5/2.54)
pie(proportions_Scer_pangenome, 
    labels = names(proportions_Scer_pangenome), 
    main = expression(italic("S. cerevisiae") ~ "'s pangenome structure"),
    col = c("red3","dodgerblue"))
#dev.off()
      #16 YPD yeast species
proportions_16_YPD_species_pangenome <- c(proportion_accessory_genes_in_the_21_species_pangenome,proportion_core_genes_in_the_21_species_pangenome)
names(proportions_16_YPD_species_pangenome) <- c("Accessory\nGenome", "Core\nGenome")
#svg(filename = paste0(output_workspace,"Yeast_subphylum_pangenome_structure.svg"),width = 7.5/2.54,height = 7.5/2.54)
pie(proportions_16_YPD_species_pangenome, 
    labels = names(proportions_16_YPD_species_pangenome), 
    main = "Yeast subphylum's pangenome\nstructure (n = 16 species)",
    col = c("red3","dodgerblue"))
#dev.off()    

#Function and orthogroup size for accessory orthogroups that are in clusters of interest (increase in expression level in time, species-specific or low expression levels)
  #Orthogroups with increase in expression level in time
v_cor_expression_level_vs_phylo_time <- rep(NA, sum(v_extract_non_problematic_orthogroups))
v_pval_expression_level_vs_phylo_time <- rep(NA, sum(v_extract_non_problematic_orthogroups))
v_proximity_to_root <- 1/tip_distances[gsub(pattern = " ",replacement = "_",x = rownames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups]),fixed = T)]
for (j in 1:sum(v_extract_non_problematic_orthogroups)){
  v_current_orthogroup_expression <- mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups][,j]
  lm_current_orthogroup_expr_vs_phylo_Time <- lm(v_current_orthogroup_expression~v_proximity_to_root)
  v_cor_expression_level_vs_phylo_time[j] <- cor(v_current_orthogroup_expression,v_proximity_to_root)
  v_pval_expression_level_vs_phylo_time[j]<- coef(summary(lm_current_orthogroup_expr_vs_phylo_Time))["v_proximity_to_root","Pr(>|t|)"]
}
v_fdr_expression_level_vs_phylo_time <- p.adjust(v_pval_expression_level_vs_phylo_time,"fdr")
v_lst_orthogroups_increasing_in_expr_with_time <- colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups])[(v_fdr_expression_level_vs_phylo_time<0.05)&(v_cor_expression_level_vs_phylo_time>0.1)]
v_lst_orthogroups_not_increasing_in_expr_with_time <- colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups])[!((v_fdr_expression_level_vs_phylo_time<0.05)&(v_cor_expression_level_vs_phylo_time>0.1))]
      #R^2 of this relationship for the orthogroups following it with significance
summary(v_cor_expression_level_vs_phylo_time[(v_fdr_expression_level_vs_phylo_time<0.05)&(v_cor_expression_level_vs_phylo_time>0.1)]^2)
      #test whether or not a pathway is enriched among orthogroups_increasing_in_expr_with_time
v_pvals_pathway_enriched_in_orthogroups_increasing_in_expr_with_time <- rep(NA, length(v_lst_KO_pathways_YPD))
names(v_pvals_pathway_enriched_in_orthogroups_increasing_in_expr_with_time) <- v_lst_KO_pathways_YPD
for (id_current_ko_pathway in 1:length(v_lst_KO_pathways_YPD)){
  current_ko_pathway <- v_lst_KO_pathways_YPD[id_current_ko_pathway]
  v_current_pathway_orthogroups <- sort(unique(subset(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY,ko_pathway==current_ko_pathway)$Orthogroup))
  mtx_cont_current_pathway <- matrix(NA,2,2)
  rownames(mtx_cont_current_pathway) <- c(paste0("Not ",current_ko_pathway),current_ko_pathway)
  colnames(mtx_cont_current_pathway) <- c("No significant increase","Significant increase")
  mtx_cont_current_pathway[1,1] <- length(intersect(v_lst_orthogroups_increasing_in_expr_with_time,v_current_pathway_orthogroups))
  mtx_cont_current_pathway[2,1] <- length(setdiff(v_lst_orthogroups_increasing_in_expr_with_time,v_current_pathway_orthogroups))
  mtx_cont_current_pathway[1,2] <- length(intersect(v_lst_orthogroups_not_increasing_in_expr_with_time,v_current_pathway_orthogroups))
  mtx_cont_current_pathway[2,2] <- length(setdiff(v_lst_orthogroups_not_increasing_in_expr_with_time,v_current_pathway_orthogroups))
  v_pvals_pathway_enriched_in_orthogroups_increasing_in_expr_with_time[id_current_ko_pathway] <- fisher.test(mtx_cont_current_pathway,alternative = "greater")$p.value
}
v_adj_pvals_fdr_pathway_enriched_in_orthogroups_increasing_in_expr_with_time <- p.adjust(v_pvals_pathway_enriched_in_orthogroups_increasing_in_expr_with_time,"fdr")
v_adj_pvals_fdr_pathway_enriched_in_orthogroups_increasing_in_expr_with_time
v_lst_ko_pathways_enriched_in_orthogroups_increasing_in_expr_with_time <- v_KO_pathway_to_pathway_title[names(v_adj_pvals_fdr_pathway_enriched_in_orthogroups_increasing_in_expr_with_time[v_adj_pvals_fdr_pathway_enriched_in_orthogroups_increasing_in_expr_with_time<0.05])]
v_lst_ko_pathways_enriched_in_orthogroups_increasing_in_expr_with_time

    #Pathway enrichment analysis for Orthogroups that are highly expressed in a specific clade, in multiple clades or both
the_mtx <- mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[,v_extract_non_problematic_orthogroups]
global_mean <- mean(the_mtx, na.rm = TRUE)
n_species <- nrow(the_mtx)
pct_above_mean <- colMeans(the_mtx > global_mean, na.rm = TRUE)
pct_is_zero <- colMeans(the_mtx == 0, na.rm = TRUE)
filtered_Orthogroups <- colnames(the_mtx)[
  pct_above_mean >= 0.10 & 
    pct_above_mean <= 0.50 & 
    pct_is_zero >= 0.50
]
#names(clade_vector) <- gsub("_", " ", names(v_clade_color[v_order_almost_21_spcs_YPD[rownames(as.matrix(the_mtx))]]), fixed = TRUE)
row_clades <- v_order_almost_21_spcs_YPD[rownames(the_mtx)]
sub_mtx <- the_mtx[, filtered_Orthogroups, drop = FALSE]
Orthogroup_expressed_in_clade <- aggregate(sub_mtx > 0, by = list(Clade = row_clades), FUN = any, na.rm = TRUE)
clade_counts_per_Orthogroup <- colSums(Orthogroup_expressed_in_clade[, -1, drop = FALSE])
clade_specific_Orthogroups <- names(clade_counts_per_Orthogroup)[clade_counts_per_Orthogroup == 1]
multi_clade_Orthogroups <- names(clade_counts_per_Orthogroup)[clade_counts_per_Orthogroup > 1]
all_Orthogroups <- colnames(the_mtx)

find_enriched_pathways <- function(query_Orthogroups, your_df, background_Orthogroups) {
  unique_pathways <- sort(unique(your_df$ko_pathway))
  the_pvals <- c()
  for (path in unique_pathways) {
    path_Orthogroups <- your_df$Orthogroup[your_df$ko_pathway == path]
    
    a <- sum(query_Orthogroups %in% path_Orthogroups)
    b <- length(query_Orthogroups) - a
    c <- sum((!(background_Orthogroups %in% query_Orthogroups)) & (background_Orthogroups %in% path_Orthogroups))
    d <- sum((!(background_Orthogroups %in% query_Orthogroups)) & (!(background_Orthogroups %in% path_Orthogroups)))
    
    contingency_matrix <- matrix(c(a, c, b, d), nrow = 2)
    p_val <- fisher.test(contingency_matrix, alternative = "greater")$p.value
    the_pvals <- c(the_pvals,p_val)
  }
  the_pvals <- p.adjust(the_pvals,"fdr")
  #names(the_pvals) <- unique_pathways
  #print(the_pvals)
  #names(the_pvals) <- NULL
  enriched_pathways <- c()
  for (i in 1:length(unique_pathways)) {
    if (!is.na(the_pvals[i]) & (the_pvals[i] < 0.05)) {
      enriched_pathways <- c(enriched_pathways, path)
    }
  }
  
  return(enriched_pathways)
}
the_df <- unique(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY[,c("Orthogroup","ko_pathway","pathway_title")])
v_pathways_enriched_clade_specific <- find_enriched_pathways(clade_specific_Orthogroups, the_df, all_Orthogroups)
sort(table( subset(the_df,Orthogroup %in% clade_specific_Orthogroups)$pathway_title ))
any(v_pathways_enriched_combined<0.05)
v_pathways_enriched_multi_clade <- find_enriched_pathways(multi_clade_Orthogroups, the_df, all_Orthogroups)
sort(table( subset(the_df,Orthogroup %in% multi_clade_Orthogroups)$pathway_title ))
any(v_pathways_enriched_multi_clade<0.05)
v_pathways_enriched_combined <- find_enriched_pathways(filtered_Orthogroups, the_df, all_Orthogroups)
sort(table( subset(the_df,Orthogroup %in% all_Orthogroups)$pathway_title ))
any(v_pathways_enriched_clade_specific<0.05)

  ###Lists of pathways evolving under positive selection vs diversifying slection vs negative selection
    #Positive selection
subset(df_ko_pathways_evo_force,pval_evo_force_call<0.05&most_likely_evo_force=="Positive selection")$pathway_title
    #Diversifying selection
subset(df_ko_pathways_evo_force,pval_evo_force_call<0.05&most_likely_evo_force=="Diversifying selection")$pathway_title
    #Negative selection
subset(df_ko_pathways_evo_force,pval_evo_force_call<0.05&most_likely_evo_force=="Negative selection")$pathway_title
    #0 pathway evolve neutrally at the transcriptomic level

###########################characterize the evolution of  paralogs ############################################################################
#Genes that have paralogs are listed in df_per_species_orthgroups_gene_features 
#Expression data are in df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY

  #Create the matrices required to analyze the evolution of paralogs at the transcriptomic level
  #control expression for differences in species transcriptome size
  #standardize fpkm at column (species) level because we are not using the R2 anymore but the expression data directly and each species not only has different library size (which is controlled for ) but also different total expression
df_yeasts_pantranscriptome$gene_X_orthogroup <- paste0(df_yeasts_pantranscriptome$Geneid,"_X_",df_yeasts_pantranscriptome$Orthogroup)
mtx_species_gene_expr_FPKM_YPD <- reshape2::acast(subset(df_yeasts_pantranscriptome,condition=="YPD"), species_smpl_lbl ~ gene_X_orthogroup, value.var = "FPKM",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)
mtx_Standardized_species_genes_FPKM_YPD <- t(scale(t(mtx_species_gene_expr_FPKM_YPD), center = TRUE, scale = TRUE))
df_Standardized_species_genes_FPKM_YPD <- reshape2::melt(mtx_Standardized_species_genes_FPKM_YPD)
colnames(df_Standardized_species_genes_FPKM_YPD) <- c("species_smpl_lbl","gene_X_orthogroup","standardized_FPKM")
df_Standardized_species_genes_FPKM_YPD$Orthogroup <- stringr::str_extract(df_Standardized_species_genes_FPKM_YPD$gene_X_orthogroup, "(?<=_X_).+")
    #mean orthogroup standardized expression 
mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD <- log2(reshape2::acast(df_Standardized_species_genes_FPKM_YPD, species_smpl_lbl ~ Orthogroup, value.var = "standardized_FPKM",fun.aggregate = function(x) mean(x,na.rm=T),fill = 0)+1)
    #median orthogroup standardized expression 
mtx_species_Orthogroups_log2_1p_median_standardized_FPKM_YPD <- log2(reshape2::acast(df_Standardized_species_genes_FPKM_YPD, species_smpl_lbl ~ Orthogroup, value.var = "standardized_FPKM",fun.aggregate = function(x) median(x,na.rm=T),fill = 0)+1)
    #sum orthogroup standardized expression 
mtx_species_Orthogroups_log2_p_cst_sum_standardized_FPKM_YPD <- log2(reshape2::acast(df_Standardized_species_genes_FPKM_YPD, species_smpl_lbl ~ Orthogroup, value.var = "standardized_FPKM",fun.aggregate = function(x) sum(x,na.rm=T),fill = 0)+abs(min(reshape2::acast(df_Standardized_species_genes_FPKM_YPD, species_smpl_lbl ~ Orthogroup, value.var = "standardized_FPKM",fun.aggregate = function(x) sum(x,na.rm=T),fill = 0),na.rm=T))+1)
    #orthogroup size
df_orthogroup_size_species_seq_in_YPD <- subset(df_yeasts_pantranscriptome,condition=="YPD") %>% group_by(species_smpl_lbl,Orthogroup) %>% summarise(orthogroup_size=length(unique(Geneid)))
mtx_orthogroup_size_in_species_seq_in_YPD <- reshape2::acast(df_orthogroup_size_species_seq_in_YPD, species_smpl_lbl ~ Orthogroup, value.var = "orthogroup_size",fun.aggregate = function(x) x,fill = 0)
mtx_dist_orthogroup_size_species_seq_in_YPD <- as.matrix(vegan::vegdist(mtx_orthogroup_size_in_species_seq_in_YPD,method = "bray"))

    #16 Yeasts species Gene content distance vs Paralogs distance (Orthogroup size comparison and ignoring double zeros) : Is duplication a major driver of gene content change
mtx_dist_species_pres_abs_YPD <- as.matrix(vegan::vegdist(mtx_gene_INTEGER_pres_abs_in_species_genomes_YPD[rownames(mtx_pairwise_expression_similarity_YPD),],method = "gower"))
      #Define the orthogroup size increase (0 if same or smaller size)
mtx_orthogroup_size_increase_for_species_YPD <- matrix(0,nrow=nrow(mtx_orthogroup_size_in_species_seq_in_YPD),ncol=nrow(mtx_orthogroup_size_in_species_seq_in_YPD))
rownames(mtx_orthogroup_size_increase_for_species_YPD) <- rownames(mtx_orthogroup_size_in_species_seq_in_YPD)
colnames(mtx_orthogroup_size_increase_for_species_YPD) <- rownames(mtx_orthogroup_size_in_species_seq_in_YPD) #this is correct because it's a square matrix
v_lst_all_Orthogroups_YPD <- unique(df_orthogroup_size_species_seq_in_YPD$Orthogroup)
nb_times_orthogroup_increased <- 0
nb_times_orthogroup_decreased <- 0
nb_times_orthogroup_stagnated <- 0
v_increase_orthogroup_size <- NULL
v_decrease_orthogroup_size <- NULL
for (i in 1:nrow(mtx_orthogroup_size_increase_for_species_YPD)){
  current_sp1 <- rownames(mtx_orthogroup_size_increase_for_species_YPD)[i]
  for (j in i:ncol(mtx_orthogroup_size_increase_for_species_YPD)){
    if (i==j){
      mtx_orthogroup_size_increase_for_species_YPD[i,j] <- NA
      next()
    }
    current_sp2 <- rownames(mtx_orthogroup_size_increase_for_species_YPD)[j]
    for (current_orthogroup in v_lst_all_Orthogroups_YPD){
      is_sp1_more_ancient_than_sp2 <- unname(tip_distances[gsub(pattern = " ",replacement = "_",x = current_sp1)] < tip_distances[gsub(pattern = " ",replacement = "_",x = current_sp2)])
      diff_size_current_orthogroup <- mtx_orthogroup_size_in_species_seq_in_YPD[current_sp2,current_orthogroup] - mtx_orthogroup_size_in_species_seq_in_YPD[current_sp1,current_orthogroup]
      
      if (diff_size_current_orthogroup==0){
        nb_times_orthogroup_stagnated <- nb_times_orthogroup_stagnated + 1
      }else{
        if (is_sp1_more_ancient_than_sp2){
          if (diff_size_current_orthogroup>0){
            nb_times_orthogroup_increased <- nb_times_orthogroup_increased + 1
            v_increase_orthogroup_size <- c(v_increase_orthogroup_size,diff_size_current_orthogroup)
            mtx_orthogroup_size_increase_for_species_YPD[i,j] <- mtx_orthogroup_size_increase_for_species_YPD[i,j] + diff_size_current_orthogroup
          }else{
            nb_times_orthogroup_decreased <- nb_times_orthogroup_decreased + 1
            v_decrease_orthogroup_size <- c(v_decrease_orthogroup_size,diff_size_current_orthogroup)
            mtx_orthogroup_size_increase_for_species_YPD[i,j] <- mtx_orthogroup_size_increase_for_species_YPD[i,j] + 0
          }
          
        }else{
          if (diff_size_current_orthogroup>0){
            nb_times_orthogroup_decreased <- nb_times_orthogroup_decreased + 1
            v_decrease_orthogroup_size <- c(v_decrease_orthogroup_size,-diff_size_current_orthogroup)
            mtx_orthogroup_size_increase_for_species_YPD[i,j] <- mtx_orthogroup_size_increase_for_species_YPD[i,j] + 0
          }else{
            nb_times_orthogroup_increased <- nb_times_orthogroup_increased + 1
            v_increase_orthogroup_size <- c(v_increase_orthogroup_size,-diff_size_current_orthogroup)
            mtx_orthogroup_size_increase_for_species_YPD[i,j] <- mtx_orthogroup_size_increase_for_species_YPD[i,j] + diff_size_current_orthogroup
          }
        }
      }
    }
  }
}
      #validations
all((v_decrease_orthogroup_size<0 )& (v_increase_orthogroup_size>0))
      #linear model gene content distance vs orthogroup set size increase (decrease and stagnation = increase of 0)
lm_YPD_species_gene_content_distance_vs_increase_orthogroup_size <- lm(get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_dist_species_pres_abs_YPD)~get_v_uppr_tri_no_diag_from_mtx(the_mtx = mtx_orthogroup_size_increase_for_species_YPD))
summary(lm_YPD_species_gene_content_distance_vs_increase_orthogroup_size)

      #compare the absolute magnitude of increases to decreases of orthogroup size
wilcox.test(abs(v_increase_orthogroup_size),abs(v_decrease_orthogroup_size))
summary(v_increase_orthogroup_size)
summary(v_decrease_orthogroup_size)
      #compare freq decrease, stagnate, increase
print(c(nb_times_orthogroup_decreased,nb_times_orthogroup_stagnated,nb_times_orthogroup_increased))

##############Develop paralogs-related metrics that are markers of orthogroups evolution################

  #Orthogroup across-species set size matrix based on df_yeasts_pantranscriptome
mtx_species_Orthogroups_set_size_YPD <- reshape2::acast(df_yeasts_pantranscriptome, species_smpl_lbl ~ Orthogroup, value.var = "Geneid",fun.aggregate = function(x) length(unique(x)),fill = 0)
#mtx_species_Orthogroups_set_size_YPD[mtx_species_Orthogroups_set_size_YPD==0] <- NA

  #Create a matrix of time progression (increase in phylogenetic Matrix of Species distance to root differences
mtx_diff_dist_to_root_for_YPD_spcs_pairs <- matrix(NA,nrow = nrow(mtx_phylo_dist_the_21_sps_YPD),ncol = ncol(mtx_phylo_dist_the_21_sps_YPD))
rownames(mtx_diff_dist_to_root_for_YPD_spcs_pairs) <- rownames(mtx_phylo_dist_the_21_sps_YPD) #gsub(pattern = "_",replacement = " ",rownames(mtx_phylo_dist_the_21_sps_YPD))
colnames(mtx_diff_dist_to_root_for_YPD_spcs_pairs) <- colnames(mtx_phylo_dist_the_21_sps_YPD) #gsub(pattern = "_",replacement = " ",colnames(mtx_phylo_dist_the_21_sps_YPD))
for (i in 1:nrow(mtx_phylo_dist_the_21_sps_YPD)){
  the_current_sp1 <- rownames(mtx_diff_dist_to_root_for_YPD_spcs_pairs)[i]
  for (j in 1:ncol(mtx_phylo_dist_the_21_sps_YPD)){
    the_current_sp2 <- colnames(mtx_diff_dist_to_root_for_YPD_spcs_pairs)[j]
    mtx_diff_dist_to_root_for_YPD_spcs_pairs[i,j] <- tip_distances[the_current_sp2] - tip_distances[the_current_sp1]
  }
}
mtx_time_progression_between_YPD_spcs_pairs <- sign(mtx_diff_dist_to_root_for_YPD_spcs_pairs)

#****PERFORM CORRELATIONS WITH MAXIMUM TIME PROGRESSION (INCREASE IN PHYLOGENETIC DISTANCE BETWEEN COLUMN Sp (sp2) and ROW Sp (Sp1)) AND FILTER OUT NEGATIVE TIME PROGRESSIONS

#Create a dataframe with the Orthogroups, the detected paralogs_fate, the power of linear regression (based on the number of species pairs sharing the orthogroup, which has been calculated above) and power bins

#Can you compute transcriptome correlation across species for a particular orthogroup

#Summarize a vectorized version of pairwise truncation scores 

#Decide of a threshold that define a set enriched in pseudogenes

#Create a dataframe with the pseudogenized Orthogroups and the detected 


#plot Pairwise difference in log2(mean_std_expr) VS Truncation score

#plot Pairwise difference in log2(mean_std_expr) VS Truncation score

##########################################Extra summary plot ideas from lab meeting
###############################################################################
#Barplots: Top 5% and Bottom 5% most/least expressed orthogroups coloured by (1) KEGG pathway title and (2) KEGG pathway class
#
#Generate color palette
cb_palette_n <- function(n) {
  base <- c(
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
    "#DDCC77", "#CC6677", "#882255", "#AA4499", "#EE8866",
    "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988"
  )
  if (n == 0) return(character(0))
  if (n <= length(base)) return(base[seq_len(n)])
  colorRampPalette(base)(n)
}

#Median expression per orthogroup across species 
v_median_expr <- apply(
  mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD,
  MARGIN = 2,        #columns = orthogroups
  FUN    = function(x) median(x, na.rm = TRUE)
)

#Top and bottom 5% orthogroups
n_pct <- ceiling(length(v_median_expr) * 0.05)
v_sorted_decr <- sort(v_median_expr, decreasing = TRUE)
df_top5 <- data.frame(
  Orthogroup        = names(v_sorted_decr)[seq_len(n_pct)],
  median_expression = unname(v_sorted_decr)[seq_len(n_pct)],
  stringsAsFactors  = FALSE
)
df_bottom5 <- data.frame(
  Orthogroup        = names(v_sorted_decr)[seq(length(v_sorted_decr),
                                               length(v_sorted_decr) - n_pct + 1)],
  median_expression = unname(v_sorted_decr)[seq(length(v_sorted_decr),
                                                length(v_sorted_decr) - n_pct + 1)],
  stringsAsFactors  = FALSE
)
#df_bottom5 is already in decreasing order of median expression

#KEGG pathway annotation per orthogroup
df_og_kegg_raw <- unique(
  df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY[
    !is.na(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY$ko_pathway),
    c("Orthogroup", "ko_pathway", "pathway_title")
  ]
)

#If pathway_class is already a column in the dataframe, use it directly;
#otherwise derive it from the ko_pathway numeric prefix (standard KEGG logic).
if ("pathway_class" %in% colnames(
  df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY)) {
  
  df_og_kegg_raw <- unique(
    df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY[
      !is.na(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY$ko_pathway),
      c("Orthogroup", "ko_pathway", "pathway_title", "pathway_class")
    ]
  )
  
}

#Primary pathway per orthogroup (first occurrence when multiple pathways exist)
df_og_kegg <- df_og_kegg_raw %>%
  group_by(Orthogroup) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    pathway_label = ifelse(!is.na(pathway_title) & nchar(pathway_title) > 0,
                           pathway_title, ko_pathway)
  )

#Join expression + KEGG data
join_kegg <- function(df) {
  out <- left_join(df, df_og_kegg, by = "Orthogroup")
  out$pathway_label[is.na(out$pathway_label)] <- "Unannotated"
  out$pathway_class[is.na(out$pathway_class)] <- "Unannotated"
  out
}

df_top5    <- join_kegg(df_top5)
df_bottom5 <- join_kegg(df_bottom5)

#Enforce x-axis factor order (decreasing median expression)
df_top5$Orthogroup <- factor(
  df_top5$Orthogroup,
  levels = df_top5$Orthogroup[order(df_top5$median_expression, decreasing = TRUE)]
)

df_bottom5$Orthogroup <- factor(
  df_bottom5$Orthogroup,
  levels = df_bottom5$Orthogroup[order(df_bottom5$median_expression, decreasing = TRUE)]
)

#Shared ggplot theme
#x-axis label size scales with number of bars
barplot_theme <- function(n_bars) {
  lbl_size <- dplyr::case_when(
    n_bars > 400 ~ 2.2,
    n_bars > 200 ~ 3.0,
    n_bars > 100 ~ 4.0,
    TRUE         ~ 5.5
  )
  theme_classic(base_size = 11) +
    theme(
      axis.text.x        = element_text(angle = 90, hjust = 1,
                                        vjust = 0.5, size = lbl_size),
      axis.title         = element_text(size = 11),
      plot.title         = element_text(size = 12, face = "bold"),
      plot.subtitle      = element_text(size = 8.5, colour = "grey40"),
      legend.title       = element_text(size = 9, face = "bold"),
      legend.text        = element_text(size = 7),
      legend.key.size    = unit(0.32, "cm"),
      legend.position    = "right",
      panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.35)
    )
}

#Generic barplot builder
make_barplot <- function(df, fill_var, fill_label,
                         title, subtitle,
                         y_label = expression(
                           "Median log"[2]*"(1 + standardized FPKM) across species")) {
  lvls <- sort(unique(df[[fill_var]]))
  pal  <- cb_palette_n(length(lvls))
  names(pal) <- lvls
  
  ggplot(df, aes(x = Orthogroup,
                 y = median_expression,
                 fill = .data[[fill_var]])) +
    geom_col(colour = NA, width = 0.85) +
    scale_fill_manual(
      values = pal,
      name   = fill_label,
      guide  = guide_legend(ncol = 1, byrow = TRUE,
                            override.aes = list(size = 3))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(title    = title,
         subtitle = subtitle,
         x        = "Orthogroup",
         y        = y_label) +
    barplot_theme(nrow(df))
}

#Build the 4 plots
lbl_top    <- paste0("n = ", nrow(df_top5),    " orthogroups")
lbl_bottom <- paste0("n = ", nrow(df_bottom5), " orthogroups")

p1 <- make_barplot(
  df          = df_top5,
  fill_var    = "pathway_label",
  fill_label  = "KEGG pathway",
  title       = "Top 5% most highly expressed orthogroups",
  subtitle    = paste(lbl_top, "| fill = KEGG pathway")
)

p2 <- make_barplot(
  df          = df_top5,
  fill_var    = "pathway_class",
  fill_label  = "KEGG pathway class",
  title       = "Top 5% most highly expressed orthogroups",
  subtitle    = paste(lbl_top, "| fill = KEGG pathway class")
)

p3 <- make_barplot(
  df          = df_bottom5,
  fill_var    = "pathway_label",
  fill_label  = "KEGG pathway",
  title       = "Bottom 5% least expressed orthogroups",
  subtitle    = paste(lbl_bottom, "| fill = KEGG pathway")
)

p4 <- make_barplot(
  df          = df_bottom5,
  fill_var    = "pathway_class",
  fill_label  = "KEGG pathway class",
  title       = "Bottom 5% least expressed orthogroups",
  subtitle    = paste(lbl_bottom, "| fill = KEGG pathway class")
)

#BOXPLOT: median expression per orthogroup by KEGG pathway CLASS
#Compute per-class median to determine x-axis order (decreasing)

df_expr_long <- df_og_kegg_raw %>%
  mutate(
    median_expr    = v_median_expr[Orthogroup],
    pathway_label  = ifelse(!is.na(pathway_title) & nchar(pathway_title) > 0,
                            pathway_title, ko_pathway),
    pathway_class  = ifelse(is.na(pathway_class), "Unannotated", pathway_class)
  ) %>%
  filter(!is.na(median_expr))

df_class_order <- df_expr_long %>%
  group_by(pathway_class) %>%
  summarise(class_median = median(median_expr, na.rm = TRUE),
            n_og = n_distinct(Orthogroup),
            .groups = "drop") %>%
  arrange(desc(class_median))

df_expr_long$pathway_class <- factor(
  df_expr_long$pathway_class,
  levels = df_class_order$pathway_class   #decreasing median left -> right
)

#Colorblind-friendly palette (Paul Tol muted, one colour per class)
n_classes  <- nrow(df_class_order)
tol_muted  <- c("#332288","#88CCEE","#44AA99","#117733","#999933",
                "#DDCC77","#CC6677","#882255","#AA4499","#EE8866")
pal_class  <- setNames(
  if (n_classes <= length(tol_muted)) tol_muted[seq_len(n_classes)]
  else colorRampPalette(tol_muted)(n_classes),
  df_class_order$pathway_class
)

#x-axis label: "Class name\n(n = X OGs)"
df_expr_long <- df_expr_long %>%
  left_join(df_class_order[, c("pathway_class","n_og")], by = "pathway_class") %>%
  mutate(class_label = paste0(pathway_class, "\n(n = ", n_og, " OGs)"))

#Factor order must match
lbl_order <- df_class_order %>%
  left_join(
    df_expr_long %>% distinct(pathway_class, n_og, class_label),
    by = c("pathway_class","n_og")
  ) %>%
  pull(class_label)

df_expr_long$class_label <- factor(df_expr_long$class_label, levels = lbl_order)

ggplot(df_expr_long,
                  aes(x = class_label,
                      y = median_expr,
                      fill = pathway_class)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4,
               linewidth = 0.45, width = 0.65) +
  stat_compare_means(method = "kruskal",paired = F) +
  scale_fill_manual(values = pal_class, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  labs(
    title    = "Median orthogroup expression by KEGG pathway class",
    subtitle = paste0("Each point = one orthogroup's median expression across ",
                      ncol(mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD),
                      " orthogroups; ordered by decreasing class median"),
    x        = "KEGG pathway class",
    y        = expression("Median log"[2]*"(1 + standardized FPKM) across species")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x        = element_text(angle = 30, hjust = 1,
                                      vjust = 1, size = 8.5),
    axis.title         = element_text(size = 11),
    plot.title         = element_text(size = 12, face = "bold"),
    plot.subtitle      = element_text(size = 8, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.4)
  )

#df_expr_long %>% group_by(pathway_title) %>% summarise(median_of_expression_median = median(median_expr,na.rm = T)) %>% View()

ggsave(filename = "Boxplot_median_OG_expression_by_KEGG_pathway_class.png",
       path = output_workspace,
       width = 20, height = 20, units = "cm", dpi = 600)
ggsave(filename = "Boxplot_median_OG_expression_by_KEGG_pathway_class.svg",
       path = output_workspace,
       width = 20 / 2.54, height = 20 / 2.54, device = svg)

#Table median orthogroup expression by KEGG pathway
df_pathway_table <- df_expr_long %>%
  group_by(ko_pathway, pathway_label, pathway_class) %>%
  summarise(
    n_orthogroups            = n_distinct(Orthogroup),
    median_of_median_expr    = median(median_expr, na.rm = TRUE),
    mean_of_median_expr      = mean(median_expr,   na.rm = TRUE),
    sd_of_median_expr        = sd(median_expr,     na.rm = TRUE),
    min_of_median_expr       = min(median_expr,    na.rm = TRUE),
    max_of_median_expr       = max(median_expr,    na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(median_of_median_expr)) %>%   #decreasing order of median
  mutate(across(where(is.numeric), ~ round(.x, 6)))

write.table(
  x         = df_pathway_table,
  file      = paste0(output_workspace,
                     "Table_median_OG_expression_by_KEGG_pathway.tsv"),
  sep       = "\t",
  na        = "NA",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)


#Quick console preview (top 10 pathways) 
cat("\nTop 10 pathways by median orthogroup expression:\n")
print(
  df_pathway_table[seq_len(min(10, nrow(df_pathway_table))),
                   c("pathway_label", "pathway_class",
                     "n_orthogroups", "median_of_median_expr")],
  row.names = FALSE
)
###############################################################################
###############Develop paralogs-related metrics that are markers of orthogroups evolution################
#
#ASSUMES the following objects are already in the environment from the main script:
#  output_workspace  : character path ending "/"
#  mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD  rows=species cols=OGs
#  mtx_species_Orthogroups_set_size_YPD                        rows=species cols=OGs; NA=absent
#  tip_distances            named numeric vector (names = Genus_species)
#  df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY cols: Orthogroup,ko_pathway,pathway_title,...
#  df_per_species_orthgroups_gene_features cols: orthogroup,species_smpl_lbl,gene_length_in_species
#  v_lst_species_smpl_lbl_YPD  character vector of species labels (space-separated)
#
#PRODUCES (all in output_workspace/Orthogroups_new_paralogs_fate/):
#  SUPP_Table_Yeast_sps_paralogs_evo_trajectory_transcriptome.tsv
#      published format: 14 columns, one row per qualifying orthogroup
#  df_paralogs_trajectory_full_metrics.tsv
#      every intermediate metric, including the marker scoring path
#  SUPP_Table_setsize_threshold_sensitivity.tsv
#  SUPP_Table_random_sign_null_by_resolution.tsv
#  SUPP_Table_candidate_marker_benchmark.tsv
#  FIG5_Expanding_orthogroups_trajectory_calls.svg / .png
#  trajectory_analysis_key_numbers.txt
#
#NOTE: the four-fate classification (Pseudogenization / Subfunctionalization /
#Neofunctionalization / Redundancy) and all of its per-fate plots have been
#removed. Orthogroups are now assigned a transcriptomic TRAJECTORY CLASS
#(Divergence / Cohesion / Unresolved) from two class-defining markers.
###############################################################################

#0. Libraries & output folders
suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(seqinr)
  library(reshape2)
  library(grid)
  library(gridExtra)
  library(ggplotify)
})

dir_fate_root <- paste0(output_workspace, "Orthogroups_new_paralogs_fate/")
dir_seq       <- paste0(output_workspace, "y1000plus_data/y1000p_orthofinder/Orthogroup_Sequences/")
if (!dir.exists(dir_fate_root)) dir.create(dir_fate_root, recursive=TRUE)

###############################################################################
#dN/dS SWITCH
#dN/dS is an INFORMATIVE METRIC ONLY; it never gates a fate call. It is
#currently DISABLED because the OrthoFinder Orthogroup_Sequences files are:
#  (1) AMINO-ACID sequences, not coding DNA  -> Ka/Ks is undefined on protein
#      (seqinr::kaks() returns a logical NA, not $ka/$ks);
#  (2) UNALIGNED and of very unequal length (e.g. 48-3504 residues);
#  (3) the ENTIRE fungal-kingdom orthogroup (~11,800 sequences across many
#      genera), not restricted to the 16 study species.
#Computing dN/dS therefore requires a SEPARATE upstream pipeline that is NOT
#part of this script:
#  - obtain the matching CDS (nucleotide) sequences for each gene,
#  - subset each orthogroup to the 16 study species,
#  - build a per-orthogroup CODON alignment (translate -> align protein with
#    MAFFT/MUSCLE -> back-translate with PAL2NAL), then
#  - run seqinr::kaks() (or codeml) on the codon alignment.
#When such codon-aligned CDS are available, set compute_dnds_enabled <- TRUE.
#While FALSE, all dN/dS columns are returned as NA, the dN/dS boxplots skip,
#and call_confidence is "supported" for every call (never "high"). NONE of the
#transcriptomic fate markers, the plurality call rule, or the truncation
#metrics are affected by this switch.
compute_dnds_enabled <- FALSE
###############################################################################

#1. Helper: save plot 20x20 cm PNG 600dpi + SVG
save_plot <- function(plt, filename_noext, dir) {
  ggsave(filename = paste0(filename_noext, ".png"),
         plot = plt, path = dir,
         width = 20, height = 20, units = "cm", dpi = 600)
  ggsave(filename = paste0(filename_noext, ".svg"),
         plot = plt, path = dir,
         width = 20 / 2.54, height = 20 / 2.54, device = svg)
}

#2. Build per-species min gene-length matrix
df_orthogroup_gene_length_summary <- df_per_species_orthgroups_gene_features %>%
  group_by(orthogroup, species_smpl_lbl) %>%
  summarise(min_gene_length  = min(gene_length_in_species, na.rm = TRUE),
            .groups = "drop")

mtx_min_gene_length_YPD <- reshape2::acast(
  df_orthogroup_gene_length_summary,
  species_smpl_lbl ~ orthogroup,
  value.var    = "min_gene_length",
  fun.aggregate = function(x) mean(x, na.rm = TRUE),
  fill         = NA_real_
)

#2b. Within-species expression rank matrix (removes scale + pathway-size bias)
mtx_species_Orthogroups_RANK_YPD <- t(apply(
  mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD, 1,
  function(x) rank(x, na.last = "keep", ties.method = "average")
))
dimnames(mtx_species_Orthogroups_RANK_YPD) <-
  dimnames(mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD)

#3. KEGG pathway maps
df_og_pathway <- unique(
  df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY[
    !is.na(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY$ko_pathway),
    c("Orthogroup", "ko_pathway", "pathway_title")
  ]
)
lst_og_pathways             <- split(df_og_pathway$ko_pathway, df_og_pathway$Orthogroup)
v_all_orthogroups_with_pathway <- names(lst_og_pathways)
lst_pathway_ogs             <- split(df_og_pathway$Orthogroup, df_og_pathway$ko_pathway)

#4. Species helpers
smpl_lbl_to_tip_key <- function(sp) gsub(" ", "_", sp, fixed = TRUE)

#5. Qualifying orthogroups (>=5 species, >=3 unique set sizes)
v_all_ogs <- intersect(
  intersect(colnames(mtx_species_Orthogroups_set_size_YPD),
            colnames(mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD)),
  v_all_orthogroups_with_pathway
)
v_nb_sp_og   <- apply(mtx_species_Orthogroups_set_size_YPD[, v_all_ogs, drop=FALSE],
                      2, function(x) sum(!is.na(x) & x > 0))
v_nb_uniq_ss <- apply(mtx_species_Orthogroups_set_size_YPD[, v_all_ogs, drop=FALSE],
                      2, function(x) length(unique(x[!is.na(x) & x > 0])))
v_qualifying_ogs <- v_all_ogs[v_nb_sp_og[v_all_ogs] >= 5 & v_nb_uniq_ss[v_all_ogs] >= 3]
cat("Qualifying orthogroups:", length(v_qualifying_ogs), "\n")

#6. Species-pair table (time_progress > 0 only)
v_sps_YPD  <- v_lst_species_smpl_lbl_YPD
v_sps_keys <- smpl_lbl_to_tip_key(v_sps_YPD)

df_all_pairs <- do.call(rbind, lapply(combn(seq_along(v_sps_YPD), 2, simplify=FALSE), function(idx) {
  sp1 <- v_sps_YPD[idx[1]]; sp2 <- v_sps_YPD[idx[2]]
  k1  <- v_sps_keys[idx[1]]; k2  <- v_sps_keys[idx[2]]
  d1  <- tip_distances[k1];   d2  <- tip_distances[k2]
  if (is.na(d1) || is.na(d2) || d1 == d2) return(NULL)
  if (d1 < d2) data.frame(sp_ancient=sp1, sp_derived=sp2,
                          dist_root_ancient=d1, dist_root_derived=d2,
                          time_progress=d2-d1, stringsAsFactors=FALSE)
  else          data.frame(sp_ancient=sp2, sp_derived=sp1,
                           dist_root_ancient=d2, dist_root_derived=d1,
                           time_progress=d1-d2, stringsAsFactors=FALSE)
}))
df_all_pairs <- df_all_pairs[df_all_pairs$time_progress > 0, ]

#7. R^2 with same / different pathway orthogroups
compute_og_R2_with_pathway <- function(focal_og, expr_mtx) {
  focal_pathways <- lst_og_pathways[[focal_og]]
  if (is.null(focal_pathways) || length(focal_pathways) == 0) return(NULL)
  v_focal <- expr_mtx[, focal_og]
  v_same  <- intersect(setdiff(unique(unlist(lst_pathway_ogs[focal_pathways])), focal_og),
                       colnames(expr_mtx))
  v_diff  <- intersect(setdiff(colnames(expr_mtx), c(focal_og, v_same)),
                       v_all_orthogroups_with_pathway)
  get_R2s <- function(ogs) {
    if (!length(ogs)) return(numeric(0))
    sapply(ogs, function(og) {
      y <- expr_mtx[, og]; vld <- !is.na(v_focal) & !is.na(y)
      if (sum(vld) < 4) NA_real_ else cor(v_focal[vld], y[vld])^2
    })
  }
  R2s  <- get_R2s(v_same); R2d <- get_R2s(v_diff)
  list(avg_R2_same = mean(R2s, na.rm=TRUE), avg_R2_diff = mean(R2d, na.rm=TRUE),
       max_R2_same = if (length(R2s)) max(R2s, na.rm=TRUE) else NA_real_,
       max_R2_diff = if (length(R2d)) max(R2d, na.rm=TRUE) else NA_real_)
}

#8. Sliding R2 / R2-ratio slope vs min set-size threshold
#Per-threshold species minimum relaxed from 4 to 2 (MIN_SP_PER_THRESHOLD) so far
#more orthogroups yield a computable slope. The user's global filters
#(>=5 species present, >=3 unique set sizes) are unchanged. If fewer than 3
#distinct set-size thresholds qualify, FALL BACK to a single Spearman of a
#per-species pathway-fidelity score against that species' raw set size, so the
#sliding markers do not become the binding NA constraint.
MIN_SP_PER_THRESHOLD <- 2

sliding_R2_slope <- function(focal_og, expr_mtx, setsize_mtx, ratio = FALSE) {
  fps   <- lst_og_pathways[[focal_og]]; if (is.null(fps)) return(list(rho=NA, pval=NA))
  v_same <- intersect(setdiff(unique(unlist(lst_pathway_ogs[fps])), focal_og), colnames(expr_mtx))
  v_diff <- intersect(setdiff(intersect(colnames(expr_mtx), v_all_orthogroups_with_pathway),
                              c(focal_og, v_same)), colnames(expr_mtx))
  if (!length(v_same)) return(list(rho=NA, pval=NA))
  
  metric_at_T <- function(T) {
    sps <- rownames(setsize_mtx)[!is.na(setsize_mtx[, focal_og]) & setsize_mtx[, focal_og] >= T]
    if (length(sps) < MIN_SP_PER_THRESHOLD) return(NA_real_)
    fe <- expr_mtx[sps, focal_og]
    R2s <- sapply(v_same, function(og) {
      y <- expr_mtx[sps, og]; vld <- !is.na(fe)&!is.na(y)
      if (sum(vld) < MIN_SP_PER_THRESHOLD) NA_real_ else cor(fe[vld],y[vld])^2 })
    avg_s <- mean(R2s, na.rm=TRUE)
    if (!ratio) return(avg_s)
    R2d <- sapply(v_diff[seq_len(min(length(v_diff),200))], function(og) {
      y <- expr_mtx[sps,og]; vld <- !is.na(fe)&!is.na(y)
      if (sum(vld) < MIN_SP_PER_THRESHOLD) NA_real_ else cor(fe[vld],y[vld])^2 })
    avg_d <- mean(R2d, na.rm=TRUE)
    if (is.na(avg_d) || avg_d==0) NA_real_ else avg_s/avg_d
  }
  
  v_ss  <- sort(unique(na.omit(setsize_mtx[, focal_og]))); v_ss <- v_ss[v_ss >= 1]
  thr   <- v_ss[-length(v_ss)]
  v_met <- if (length(thr)) sapply(thr, metric_at_T) else numeric(0)
  vld   <- !is.na(v_met)
  
  #Record WHICH statistic produced the returned rho.
  #thr = v_ss[-length(v_ss)] drops the highest set size, so the number of
  #available thresholds is (number of distinct NON-ZERO set sizes) - 1. The
  #guard below therefore requires >= 4 distinct non-zero set sizes before the
  #primary sliding-threshold statistic can be used. Below that the function
  #silently switches to a DIFFERENT statistic (the per-species fallback), which
  #is not a degraded version of the same test. The path is now returned so that
  #downstream calls can be stratified by it instead of being reconstructed.
  if (sum(vld) >= 3) {
    ct <- cor.test(thr[vld], v_met[vld], method="spearman", exact=FALSE)
    return(list(rho=unname(ct$estimate), pval=ct$p.value,
                path="threshold", n_thr=length(thr), n_valid=sum(vld)))
  }
  
  #Fallback: per-species pathway-fidelity vs that species' raw set size
  sps_all <- rownames(setsize_mtx)[!is.na(setsize_mtx[, focal_og]) & setsize_mtx[, focal_og] >= 1]
  if (length(sps_all) < 4) return(list(rho=NA, pval=NA))
  ss_vec    <- setsize_mtx[sps_all, focal_og]
  same_expr <- expr_mtx[sps_all, v_same, drop=FALSE]
  centroid  <- rowMeans(same_expr, na.rm=TRUE)
  fe_all    <- expr_mtx[sps_all, focal_og]
  fidelity  <- -abs(fe_all - centroid)   #higher = closer to pathway centroid
  vld2 <- !is.na(fidelity) & !is.na(ss_vec)
  if (sum(vld2) < 4 || length(unique(ss_vec[vld2])) < 2) return(list(rho=NA, pval=NA))
  ct <- cor.test(ss_vec[vld2], fidelity[vld2], method="spearman", exact=FALSE)
  list(rho=unname(ct$estimate), pval=ct$p.value,
       path="fallback", n_thr=length(thr), n_valid=sum(vld))
}

#9. ESD (rank-normalised Jensen-Shannon divergence) slope
#Per-threshold species minimum relaxed to MIN_SP_PER_THRESHOLD (2). Fallback when
#fewer than 3 thresholds qualify: single Spearman of per-species JSD-to-pathway-
#centroid against that species' raw set size.
ESD_slope <- function(focal_og, rank_mtx, setsize_mtx) {
  fps  <- lst_og_pathways[[focal_og]]; if (is.null(fps)) return(list(rho=NA, pval=NA))
  same_ogs <- intersect(setdiff(unique(unlist(lst_pathway_ogs[fps])), focal_og),
                        colnames(rank_mtx))
  if (!length(same_ogs)) return(list(rho=NA, pval=NA))
  
  jsd_for_sps <- function(sps) {
    if (length(sps) < MIN_SP_PER_THRESHOLD) return(NA_real_)
    p <- rank_mtx[sps, focal_og] + 0.01
    q <- rowMeans(rank_mtx[sps, same_ogs, drop=FALSE], na.rm=TRUE) + 0.01
    p <- p/sum(p); q <- q/sum(q); m <- 0.5*(p+q)
    kl <- function(a,b) sum(a*log(a/b), na.rm=TRUE)
    sqrt(0.5*kl(p,m) + 0.5*kl(q,m))
  }
  
  v_ss <- sort(unique(na.omit(setsize_mtx[, focal_og]))); v_ss <- v_ss[v_ss >= 1]
  thr  <- v_ss[-length(v_ss)]
  v_esd <- if (length(thr)) sapply(thr, function(T) {
    sps <- rownames(setsize_mtx)[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= T]
    jsd_for_sps(sps)
  }) else numeric(0)
  vld <- !is.na(v_esd)
  
  #Sliding_R2_slope: Same threshold logic as above and same >= 4 distinct non-zero set-size requirement to reach the primary test.
  if (sum(vld) >= 3) {
    ct <- cor.test(thr[vld], v_esd[vld], method="spearman", exact=FALSE)
    return(list(rho=unname(ct$estimate), pval=ct$p.value,
                path="threshold", n_thr=length(thr), n_valid=sum(vld)))
  }
  
  #Fallback: per-species divergence-from-centroid vs raw set size
  sps_all <- rownames(setsize_mtx)[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= 1]
  if (length(sps_all) < 4) return(list(rho=NA, pval=NA))
  ss_vec   <- setsize_mtx[sps_all, focal_og]
  centroid <- rowMeans(rank_mtx[sps_all, same_ogs, drop=FALSE], na.rm=TRUE)
  diverg   <- abs(rank_mtx[sps_all, focal_og] - centroid)  #higher = more divergent
  vld2 <- !is.na(diverg) & !is.na(ss_vec)
  if (sum(vld2) < 4 || length(unique(ss_vec[vld2])) < 2) return(list(rho=NA, pval=NA))
  ct <- cor.test(ss_vec[vld2], diverg[vld2], method="spearman", exact=FALSE)
  list(rho=unname(ct$estimate), pval=ct$p.value,
       path="fallback", n_thr=length(thr), n_valid=sum(vld))
}

#10. Normalised CV slope
#Per-threshold species minimum relaxed to MIN_SP_PER_THRESHOLD (2).
#GUARD against unstable ratios: the normalised CV is only computed when the
#pathway-median CV exceeds CV_FLOOR; otherwise that threshold returns NA. The
#final ratio is winsorised to +/- CV_WINSOR to stop extreme outliers (the
#previous run produced values from -455 to +170) from dominating the slope.
#Fallback when fewer than 3 thresholds qualify: single Spearman of per-species
#normalised CV proxy is not meaningful (CV needs multiple species), so we keep
#the threshold approach and return NA if it cannot be computed.
CV_FLOOR  <- 1e-3
CV_WINSOR <- 10

norm_CV_slope <- function(focal_og, expr_mtx, setsize_mtx) {
  fps <- lst_og_pathways[[focal_og]]; if (is.null(fps)) return(list(slope=NA, pval=NA))
  same_ogs <- intersect(setdiff(unique(unlist(lst_pathway_ogs[fps])), focal_og), colnames(expr_mtx))
  if (!length(same_ogs)) return(list(slope=NA, pval=NA))
  v_ss <- sort(unique(na.omit(setsize_mtx[, focal_og]))); v_ss <- v_ss[v_ss >= 1]
  thr  <- v_ss[-length(v_ss)]
  if (!length(thr)) return(list(slope=NA, pval=NA))
  v_ncv <- sapply(thr, function(T) {
    sps <- rownames(setsize_mtx)[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= T]
    if (length(sps) < MIN_SP_PER_THRESHOLD) return(NA_real_)
    fe <- expr_mtx[sps, focal_og]
    cv_f <- sd(fe, na.rm=TRUE) / (abs(mean(fe, na.rm=TRUE)) + 1e-9)
    cvs <- sapply(same_ogs, function(og) {
      e <- expr_mtx[sps, og]; sd(e,na.rm=TRUE)/(abs(mean(e,na.rm=TRUE))+1e-9) })
    med <- median(cvs, na.rm=TRUE)
    if (is.na(med) || med < CV_FLOOR) return(NA_real_)   #guard near-zero denominator
    r <- cv_f/med
    max(min(r, CV_WINSOR), -CV_WINSOR)                    #winsorise
  })
  vld <- !is.na(v_ncv)
  if (sum(vld) < 3) return(list(slope=NA, pval=NA))
  sm <- summary(lm(v_ncv[vld] ~ thr[vld]))$coefficients
  list(slope=unname(sm[2,1]), pval=unname(sm[2,4]))
}

#11. Partial regression: expr_increase ~ setsize_increase + time-
partial_reg_expr_setsize_time <- function(focal_og, expr_mtx, setsize_mtx) {
  sps_present <- v_sps_YPD[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= 1]
  rows <- do.call(rbind, lapply(seq_len(nrow(df_all_pairs)), function(i) {
    anc <- df_all_pairs$sp_ancient[i]; der <- df_all_pairs$sp_derived[i]
    if (!(anc %in% sps_present) || !(der %in% sps_present)) return(NULL)
    ss_a <- setsize_mtx[anc, focal_og]; ss_d <- setsize_mtx[der, focal_og]
    if (is.na(ss_a)||is.na(ss_d)||ss_d-ss_a<=0) return(NULL)
    ea <- expr_mtx[anc, focal_og]; ed <- expr_mtx[der, focal_og]
    if (is.na(ea)||is.na(ed)||ed-ea<=0) return(NULL)
    data.frame(expr_increase=ed-ea, setsize_increase=ss_d-ss_a,
               time_progress=df_all_pairs$time_progress[i])
  }))
  if (is.null(rows)||nrow(rows)<5)
    return(list(beta_setsize=NA, pval_setsize=NA, beta_time=NA, pval_time=NA))
  sm <- summary(lm(expr_increase ~ setsize_increase + time_progress, data=rows))$coefficients
  list(beta_setsize = if ("setsize_increase" %in% rownames(sm)) sm["setsize_increase",1] else NA,
       pval_setsize = if ("setsize_increase" %in% rownames(sm)) sm["setsize_increase",4] else NA,
       beta_time    = if ("time_progress"    %in% rownames(sm)) sm["time_progress",   1] else NA,
       pval_time    = if ("time_progress"    %in% rownames(sm)) sm["time_progress",   4] else NA)
}

#12. Logarithmic-decay test (replaces AICc exponential comparison)
#Per user request: drop the AIC/AICc-vs-linear model comparison entirely and
#instead test whether set-size increase follows a LOGARITHMIC decay in time, i.e.
#fit  delta_ss ~ a + b * log(time_progress)  and require the log-time coefficient
#b to be significantly NEGATIVE. This is a plain linear model in log(time):
#  - always converges (no NLS failures),
#  - the slope and its p-value are read directly from the model summary,
#  - "significant negative log-time component" = decay slows as time increases.
#Returns is_log_decay (logical), the slope, and its raw p-value. Stored in the
#output as log_decay / log_decay_slope / log_decay_pval (replacing the previous
#exp_decay / nls_converged columns and the AICc model comparison entirely).
log_decay_test <- function(focal_og, setsize_mtx) {
  sps_present <- v_sps_YPD[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= 1]
  rows <- do.call(rbind, lapply(seq_len(nrow(df_all_pairs)), function(i) {
    anc <- df_all_pairs$sp_ancient[i]; der <- df_all_pairs$sp_derived[i]
    if (!(anc %in% sps_present)||(! der %in% sps_present)) return(NULL)
    ss_a <- setsize_mtx[anc, focal_og]; ss_d <- setsize_mtx[der, focal_og]
    if (is.na(ss_a)||is.na(ss_d)||ss_d-ss_a<=0) return(NULL)
    tp <- df_all_pairs$time_progress[i]
    if (is.na(tp) || tp <= 0) return(NULL)            #log(time) requires time > 0
    data.frame(delta_ss=ss_d-ss_a, time_progress=tp)
  }))
  if (is.null(rows)||nrow(rows)<5)
    return(list(is_log_decay=NA, slope=NA, pval=NA))
  fit <- tryCatch(lm(delta_ss ~ log(time_progress), data=rows),
                  error=function(e) NULL)
  if (is.null(fit)) return(list(is_log_decay=FALSE, slope=NA, pval=NA))
  sm <- summary(fit)$coefficients
  if (!("log(time_progress)" %in% rownames(sm)))
    return(list(is_log_decay=FALSE, slope=NA, pval=NA))
  b_est  <- sm["log(time_progress)", 1]
  pval_b <- sm["log(time_progress)", 4]
  is_decay <- (b_est < 0) && (pval_b < 0.05)          #significant negative log-time
  list(is_log_decay=is_decay, slope=b_est, pval=pval_b)
}

#13. Truncation metrics (universal for all 4 fates)
truncation_metrics <- function(focal_og, setsize_mtx, min_len_mtx) {
  sps_present <- v_sps_YPD[!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= 1]
  rows <- do.call(rbind, lapply(seq_len(nrow(df_all_pairs)), function(i) {
    anc <- df_all_pairs$sp_ancient[i]; der <- df_all_pairs$sp_derived[i]
    if (!(anc %in% sps_present)||(! der %in% sps_present)) return(NULL)
    ss_a <- setsize_mtx[anc, focal_og]; ss_d <- setsize_mtx[der, focal_og]
    if (is.na(ss_a)||is.na(ss_d)||ss_d<=ss_a) return(NULL)
    ml_a <- min_len_mtx[anc, focal_og]; ml_d <- min_len_mtx[der, focal_og]
    if (is.na(ml_a)||is.na(ml_d)||ml_a==0) return(NULL)
    pct <- (ml_a - ml_d) / ml_a
    data.frame(delta_ss=ss_d-ss_a, pct_decrease=pct, sign_trunc=pct>0.30)
  }))
  if (is.null(rows)||nrow(rows)==0)
    return(list(nb_pairs_sign_trunc=NA_integer_, avg_truncation_score=NA_real_))
  list(nb_pairs_sign_trunc  = sum(rows$sign_trunc, na.rm=TRUE),
       avg_truncation_score = mean(rows$pct_decrease, na.rm=TRUE))
}

#14. dN/dS: pairwise cross-species + within-species rho ~ set size
#dN/dS is purely a METRIC (not a marker for fate calling) and is DISABLED by
#default (compute_dnds_enabled <- FALSE at the top of the script) because the
#current OrthoFinder files are unaligned amino acids spanning the whole fungal
#kingdom, from which Ka/Ks cannot be computed. See the switch block above for
#the codon-aligned-CDS pipeline required to enable it.
#IMPLEMENTATION NOTE (for when enabled): seqinr::kaks() takes an `alignment`
#object and returns ka and ks as pairwise dist matrices over ALL sequences at
#once (it does NOT accept two individual sequences). We read each file with
#read.alignment(), keep only sequences mapping to a study species (strict
#genus_species match), and call kaks() ONCE on that sub-alignment.
#Returns:
#  mean_dNdS         mean pairwise dN/dS over CROSS-species pairs (Ks > 0.001)
#  max_dNdS          max  pairwise dN/dS over cross-species pairs
#  n_pairs_high_Ks   diagnostic: count of cross-species pairs with Ks > 1 (saturation)
#  rho_dNdS_setsize  partial Spearman rho between within-species mean dN/dS and
#                    set size, controlling for time_progress (>= 5 multi-copy sp.)
#
#Cross-species pairs = the two sequences come from DIFFERENT study species.
#Within-species dN/dS = mean over paralog pairs whose two sequences are from the
#SAME species (only species with set size >= 2 contribute), capturing relaxation
#of constraint as copy number grows.

compute_dNdS <- function(focal_og, setsize_mtx) {
  empty <- list(mean_dNdS=NA, max_dNdS=NA, n_pairs_high_Ks=NA,
                rho_dNdS_setsize=NA, pval_rho_dNdS=NA)
  #dN/dS is disabled unless codon-aligned CDS are available (see switch at top).
  #When disabled, return all-NA immediately without touching seqinr.
  if (!isTRUE(compute_dnds_enabled)) return(empty)
  
  fasta_path <- paste0(dir_seq, focal_og, ".fa")
  if (!file.exists(fasta_path)) fasta_path <- paste0(dir_seq, focal_og, ".fasta")
  if (!file.exists(fasta_path)) return(empty)
  
  #NOTE: kaks() requires CODON-ALIGNED NUCLEOTIDE sequences. The current
  #OrthoFinder files are unaligned amino acids spanning the whole fungal kingdom,
  #so this branch only runs meaningfully once a proper codon-aligned CDS file
  #(restricted to the 16 study species) is supplied and compute_dnds_enabled=TRUE.
  #seqinr::kaks() operates on an `alignment` object and returns ka / ks as
  #pairwise dist matrices over ALL sequences at once (it does NOT take two
  #individual sequences). We read the file as an alignment, keep only sequences
  #that map to a study species, call kaks() ONCE on that sub-alignment, and
  #index the returned matrices for the cross-species and within-species summaries.
  aln <- tryCatch(read.alignment(fasta_path, format="fasta"),
                  error=function(e) NULL)
  #read.alignment must return a list-like alignment; guard with is.list() so a
  #NULL or unexpected atomic return can never trigger "$ on atomic vector".
  if (!is.list(aln) || is.null(aln[["nb"]]) || !length(aln[["nb"]]) || aln[["nb"]] < 2)
    return(empty)
  
  #Map each alignment sequence to a study species-
  #Headers look like "geneID|Order_species.final". Match lowercase-to-lowercase
  #and require the FULL genus_species key to appear in the header. The previous
  #epithet-only fallback is REMOVED: with kingdom-wide files containing many
  #genera, a shared species epithet (e.g. "anomalus" in both
  #brettanomyces_anomalus and Wickerhamomyces_anomalus) produced false matches.
  v_sps_keys_lc  <- tolower(gsub(" ", "_", v_sps_YPD, fixed=TRUE))
  
  match_species <- function(nm) {
    nm_low <- tolower(nm)
    idx <- which(vapply(v_sps_keys_lc, function(k) grepl(k, nm_low, fixed=TRUE), logical(1)))
    if (length(idx) >= 1) return(v_sps_YPD[idx[1]])   #full genus_species only
    NA_character_
  }
  
  seq_sp_all <- vapply(aln[["nam"]], match_species, character(1))
  keep_idx   <- which(!is.na(seq_sp_all))
  if (length(keep_idx) < 2) return(empty)
  
  #Subset the alignment to matched sequences only
  sub_aln <- list(
    nb  = length(keep_idx),
    nam = aln[["nam"]][keep_idx],
    seq = aln[["seq"]][keep_idx],
    com = if (!is.null(aln[["com"]])) aln[["com"]][keep_idx] else NA
  )
  class(sub_aln) <- "alignment"
  seq_sp <- seq_sp_all[keep_idx]          #species label per retained sequence
  N <- length(keep_idx)
  
  #One kaks() call over the whole sub-alignment
  kk <- tryCatch(suppressWarnings(kaks(sub_aln)), error=function(e) NULL)
  #kaks() must return a list with ka/ks; guard with is.list() so an atomic/NULL
  #return cannot trigger "$ operator is invalid for atomic vectors".
  if (!is.list(kk) || is.null(kk[["ka"]]) || is.null(kk[["ks"]])) return(empty)
  ka_m <- tryCatch(as.matrix(kk[["ka"]]), error=function(e) NULL)
  ks_m <- tryCatch(as.matrix(kk[["ks"]]), error=function(e) NULL)
  if (is.null(ka_m) || is.null(ks_m)) return(empty)
  #Guard against degenerate matrices
  if (any(dim(ka_m) != c(N, N)) || any(dim(ks_m) != c(N, N))) return(empty)
  
  #valid pairwise (i<j) dN/dS with Ks filter (Ks > 0.001; saturation Ks>1 flagged)
  pair_dn  <- numeric(0); pair_ks <- numeric(0)
  pair_i   <- integer(0); pair_j  <- integer(0)
  for (i in seq_len(N-1)) for (j in (i+1):N) {
    ks <- ks_m[i, j]; ka <- ka_m[i, j]
    if (is.na(ks) || is.na(ka) || ks <= 0.001 || ka < 0) next
    pair_dn <- c(pair_dn, ka/ks); pair_ks <- c(pair_ks, ks)
    pair_i  <- c(pair_i, i);      pair_j  <- c(pair_j, j)
  }
  if (length(pair_dn) == 0) return(empty)
  
  #Cross-species summary (pairs whose two sequences are from DIFFERENT sp) -
  diff_sp   <- seq_sp[pair_i] != seq_sp[pair_j]
  xsp_dn    <- pair_dn[diff_sp]; xsp_ks <- pair_ks[diff_sp]
  n_high_ks <- sum(xsp_ks > 1, na.rm=TRUE)
  mean_dN   <- if (length(xsp_dn) >= 3) mean(xsp_dn) else NA_real_
  max_dN    <- if (length(xsp_dn) >= 3) max(xsp_dn)  else NA_real_
  
  #Within-species mean dN/dS per multi-copy species
  sps_multicopy <- v_sps_YPD[!is.na(setsize_mtx[, focal_og]) & setsize_mtx[, focal_og] >= 2]
  within_rows <- do.call(rbind, lapply(sps_multicopy, function(sp) {
    sel <- which(seq_sp[pair_i] == sp & seq_sp[pair_j] == sp)  #both seqs from sp
    if (!length(sel)) return(NULL)
    vals <- pair_dn[sel]
    if (!length(vals)) return(NULL)
    data.frame(species=sp,
               within_mean_dNdS = mean(vals),
               set_size         = setsize_mtx[sp, focal_og],
               time_progress    = unname(tip_distances[smpl_lbl_to_tip_key(sp)]),
               stringsAsFactors = FALSE)
  }))
  
  rho_dN <- NA_real_; pval_rho <- NA_real_
  if (!is.null(within_rows) && nrow(within_rows) >= 5 &&
      length(unique(within_rows$set_size)) >= 2) {
    #partial Spearman: rho(within_mean_dNdS, set_size | time_progress)
    rk_y <- rank(within_rows$within_mean_dNdS, ties.method="average")
    rk_x <- rank(within_rows$set_size,         ties.method="average")
    rk_t <- rank(within_rows$time_progress,    ties.method="average")
    res_y <- residuals(lm(rk_y ~ rk_t))
    res_x <- residuals(lm(rk_x ~ rk_t))
    if (sd(res_x) > 0 && sd(res_y) > 0) {
      ct <- suppressWarnings(cor.test(res_x, res_y, method="pearson"))
      rho_dN <- unname(ct$estimate); pval_rho <- ct$p.value
    }
  }
  
  list(mean_dNdS=mean_dN, max_dNdS=max_dN, n_pairs_high_Ks=n_high_ks,
       rho_dNdS_setsize=rho_dN, pval_rho_dNdS=pval_rho)
}

#15. Parallelised main loop over qualifying orthogroups
nb_cores_fate <- 6
cl <- makeCluster(nb_cores_fate)
registerDoParallel(cl)
cat("Starting parallelised orthogroup fate analysis on", nb_cores_fate, "cores...\n")

clusterExport(cl, varlist = c(
  "output_workspace","dir_seq",
  "mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD",
  "mtx_species_Orthogroups_set_size_YPD",
  "mtx_species_Orthogroups_RANK_YPD",
  "mtx_min_gene_length_YPD",
  "tip_distances","df_all_pairs","v_sps_YPD","v_sps_keys",
  "lst_og_pathways","lst_pathway_ogs","v_all_orthogroups_with_pathway",
  "smpl_lbl_to_tip_key",
  "compute_og_R2_with_pathway","sliding_R2_slope",
  "ESD_slope","norm_CV_slope","partial_reg_expr_setsize_time",
  "log_decay_test","truncation_metrics","compute_dNdS",
  "MIN_SP_PER_THRESHOLD","CV_FLOOR","CV_WINSOR"
), envir = environment())
clusterEvalQ(cl, { library(seqinr); library(dplyr) })

lst_fate_results <- foreach(
  focal_og       = v_qualifying_ogs,
  .packages      = c("seqinr"),
  .errorhandling = "pass"
) %dopar% {
  
  expr_mtx    <- mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD
  rank_mtx    <- mtx_species_Orthogroups_RANK_YPD
  setsize_mtx <- mtx_species_Orthogroups_set_size_YPD
  min_len_mtx <- mtx_min_gene_length_YPD
  
  #Safe field accessor: returns x[[field]] only when x is a list that actually
  #has it; otherwise NA. Prevents "$ operator is invalid for atomic vectors" if
  #any helper returns an atomic / NULL in an edge case.
  g <- function(x, field) {
    if (is.list(x) && !is.null(x[[field]]) && length(x[[field]])) x[[field]] else NA
  }
  
  nb_sp      <- sum(!is.na(setsize_mtx[,focal_og]) & setsize_mtx[,focal_og] >= 1)
  
  #count DISTINCT NON-ZERO set sizes only.
  #The previous version used length(unique(na.omit(setsize_mtx[,focal_og]))),
  #which counts a set size of 0 (species where the orthogroup is ABSENT) as a
  #distinct level. That inflated nb_unique_set_sizes by 1 for every orthogroup
  #missing from at least one species, and made the reported column inconsistent
  #with the qualifying filter in section 5 (which already used x > 0).
  #Consequence of the old behaviour: the reported column conflated copy-number
  #resolution with presence/absence, so accessory orthogroups appeared to have
  #one more level of resolution than they really had.
  ss_col     <- setsize_mtx[, focal_og]
  nb_uniq_ss <- length(unique(ss_col[!is.na(ss_col) & ss_col > 0]))
  
  r2_info        <- compute_og_R2_with_pathway(focal_og, expr_mtx)
  avg_R2_same    <- g(r2_info, "avg_R2_same")
  avg_R2_diff    <- g(r2_info, "avg_R2_diff")
  max_R2_same    <- g(r2_info, "max_R2_same")
  max_R2_diff    <- g(r2_info, "max_R2_diff")
  ratio_R2       <- if (!is.na(avg_R2_diff) && avg_R2_diff > 0 && !is.na(avg_R2_same))
    avg_R2_same / avg_R2_diff else NA
  
  slope_R2_same  <- sliding_R2_slope(focal_og, expr_mtx, setsize_mtx, ratio=FALSE)
  slope_R2_ratio <- sliding_R2_slope(focal_og, expr_mtx, setsize_mtx, ratio=TRUE)
  esd_sl         <- ESD_slope(focal_og, rank_mtx, setsize_mtx)
  ncv_sl         <- norm_CV_slope(focal_og, expr_mtx, setsize_mtx)
  preg           <- partial_reg_expr_setsize_time(focal_og, expr_mtx, setsize_mtx)
  log_dec        <- log_decay_test(focal_og, setsize_mtx)
  trunc          <- truncation_metrics(focal_og, setsize_mtx, min_len_mtx)
  dnds           <- compute_dNdS(focal_og, setsize_mtx)
  
  #===========================================================================
  #The loop records raw per-orthogroup metrics only. The trajectory class is
  #assigned post-hoc in section 16 from the two class-defining markers.
  #(The former four-fate plurality rule has been removed entirely.)
  #===========================================================================
  
  data.frame(
    Orthogroup                         = focal_og,
    nb_species                         = nb_sp,
    nb_unique_set_sizes                = nb_uniq_ss,
    avg_truncation_score               = g(trunc, "avg_truncation_score"),
    avg_nb_species_pairs_sign_trunc    = g(trunc, "nb_pairs_sign_trunc"),
    ESD_rho                            = g(esd_sl, "rho"),
    norm_CV_slope                      = g(ncv_sl, "slope"),
    beta_setsize                       = g(preg, "beta_setsize"),
    pval_setsize                       = g(preg, "pval_setsize"),
    beta_time                          = g(preg, "beta_time"),
    pval_time                          = g(preg, "pval_time"),
    slope_R2_same_rho                  = g(slope_R2_same, "rho"),
    slope_R2_same_pval                 = g(slope_R2_same, "pval"),
    #which statistic scored each marker, and how many thresholds existed
    R2slope_path                       = g(slope_R2_same, "path"),
    R2slope_n_thr                      = g(slope_R2_same, "n_thr"),
    ESD_path                           = g(esd_sl, "path"),
    ESD_n_thr                          = g(esd_sl, "n_thr"),
    slope_R2_ratio_rho                 = g(slope_R2_ratio, "rho"),
    slope_R2_ratio_pval                = g(slope_R2_ratio, "pval"),
    ratio_R2_same_diff                 = ratio_R2,
    max_R2_same                        = max_R2_same,
    max_R2_diff                        = max_R2_diff,
    log_decay                          = g(log_dec, "is_log_decay"),
    log_decay_slope                    = g(log_dec, "slope"),
    log_decay_pval                     = g(log_dec, "pval"),
    norm_CV_pval                       = g(ncv_sl, "pval"),
    #dN/dS metrics (informative only; not used in fate call)
    mean_dNdS                          = g(dnds, "mean_dNdS"),
    max_dNdS                           = g(dnds, "max_dNdS"),
    n_pairs_high_Ks                    = g(dnds, "n_pairs_high_Ks"),
    rho_dNdS_setsize                   = g(dnds, "rho_dNdS_setsize"),
    pval_rho_dNdS                      = g(dnds, "pval_rho_dNdS"),
    stringsAsFactors                   = FALSE
  )
}

stopCluster(cl)
registerDoSEQ()

###############################################################################
#16. Trajectory-class calls from the TWO class-defining markers
#
#The former four-fate plurality rule (Pseudogenization / Subfunctionalization /
#Neofunctionalization / Redundancy) has been REMOVED. Classes are now assigned
#from the two pathway-fidelity markers only:
#
#  ESD_rho              Spearman of the Jensen-Shannon divergence between the
#                       orthogroup's rank-expression profile and its pathway
#                       centroid, across sliding set-size thresholds.
#                       rho > 0  -> Divergence      (drifts from its pathway)
#                       rho <= 0 -> Cohesion        (stays anchored)
#
#  slope_R2_same_rho    Spearman of the mean R2 with same-pathway partners,
#                       across the same thresholds.
#                       rho < 0  -> Divergence      (loses co-regulation)
#                       rho >= 0 -> Cohesion        (retains co-regulation)
#
#A class is assigned ONLY where the two markers agree in direction. Where they
#conflict the orthogroup is left Unresolved rather than decided by a further
#marker (none of the candidates in section 20 beat a non-informative constant).
#
#TERMINOLOGY: "Cohesion" replaces the former "Convergence" throughout. Cohesion
#is the umbrella for expression states held more similar than the drift
#expectation, comprising both directional convergence and conservation.
###############################################################################
is_df_row <- vapply(lst_fate_results, is.data.frame, logical(1))
n_failed  <- sum(!is_df_row)
if (n_failed > 0) {
  cat("WARNING:", n_failed, "of", length(lst_fate_results),
      "orthogroups errored inside the worker and were dropped.\n")
  errs <- lst_fate_results[!is_df_row]
  msgs <- unique(vapply(errs, function(e) {
    if (inherits(e, "condition")) conditionMessage(e) else paste(class(e), collapse=",")
  }, character(1)))
  cat("  Distinct worker error messages (up to 5):\n")
  for (m in head(msgs, 5)) cat("   -", m, "\n")
}
if (!any(is_df_row)) stop("All orthogroups failed; no results to assemble.")
df_paralogs_fate <- do.call(rbind, lst_fate_results[is_df_row])

#marker votes
CLASS_DIV <- "Divergence"
CLASS_COH <- "Cohesion"
CLASS_UNR <- "Unresolved"

vote_esd <- function(rho) ifelse(is.na(rho), NA_character_,
                                 ifelse(rho >  0, CLASS_DIV, CLASS_COH))
vote_r2  <- function(rho) ifelse(is.na(rho), NA_character_,
                                 ifelse(rho <  0, CLASS_DIV, CLASS_COH))

df_paralogs_fate$ESD_vote     <- vote_esd(df_paralogs_fate$ESD_rho)
df_paralogs_fate$R2slope_vote <- vote_r2( df_paralogs_fate$slope_R2_same_rho)

df_paralogs_fate$markers_agree <-
  !is.na(df_paralogs_fate$ESD_vote) & !is.na(df_paralogs_fate$R2slope_vote) &
  (df_paralogs_fate$ESD_vote == df_paralogs_fate$R2slope_vote)

df_paralogs_fate$min_abs_rho <- pmin(abs(df_paralogs_fate$ESD_rho),
                                     abs(df_paralogs_fate$slope_R2_same_rho))

df_paralogs_fate$Transcriptomic_trajectory_class <- ifelse(
  df_paralogs_fate$markers_agree, df_paralogs_fate$ESD_vote, CLASS_UNR)

#why is an orthogroup unresolved?
#Detection criteria only (no results): the conflict is classed by the WEAKER of
#the two marker correlations.
#  strong_discordant_signal : both |rho| >= 0.5  -> two strong, opposing signals
#  weak_discordant_signal   : 0.2 <= min|rho| < 0.5
#  low_signal_sign_unstable : min|rho| < 0.2     -> sign is close to arbitrary
RHO_STRONG <- 0.5
RHO_WEAK   <- 0.2
df_paralogs_fate$unresolved_reason <- ifelse(
  df_paralogs_fate$markers_agree, NA_character_,
  ifelse(is.na(df_paralogs_fate$min_abs_rho), "marker_not_computable",
         ifelse(df_paralogs_fate$min_abs_rho >= RHO_STRONG, "strong_discordant_signal",
                ifelse(df_paralogs_fate$min_abs_rho >= RHO_WEAK,   "weak_discordant_signal",
                       "low_signal_sign_unstable"))))

#scoring path
#An orthogroup reaches the PRIMARY sliding-threshold statistic only when it has
#at least 3 usable thresholds, i.e. at least 4 distinct non-zero set sizes.
#Below that both markers switch to the per-species fallback statistic.
df_paralogs_fate$ESD_path     <- as.character(df_paralogs_fate$ESD_path)
df_paralogs_fate$R2slope_path <- as.character(df_paralogs_fate$R2slope_path)
df_paralogs_fate$scoring_path <- ifelse(
  df_paralogs_fate$ESD_path == "threshold" & df_paralogs_fate$R2slope_path == "threshold",
  "threshold",
  ifelse(is.na(df_paralogs_fate$ESD_path) | is.na(df_paralogs_fate$R2slope_path),
         NA_character_, "fallback"))

df_paralogs_fate$Transcriptomic_trajectory_class <- factor(
  df_paralogs_fate$Transcriptomic_trajectory_class,
  levels = c(CLASS_DIV, CLASS_COH, CLASS_UNR))

cat("\n=== Trajectory-class calls (2-marker rule) ===\n")
print(table(df_paralogs_fate$Transcriptomic_trajectory_class, useNA="ifany"))
cat("\nUnresolved reasons:\n")
print(table(df_paralogs_fate$unresolved_reason, useNA="no"))
cat("\nScoring path:\n")
print(table(df_paralogs_fate$scoring_path, useNA="ifany"))
cat("\nClass by scoring path:\n")
print(table(df_paralogs_fate$scoring_path,
            df_paralogs_fate$Transcriptomic_trajectory_class, useNA="ifany"))

#full intermediate table (all metrics retained)
write.table(df_paralogs_fate,
            file = paste0(dir_fate_root, "df_paralogs_trajectory_full_metrics.tsv"),
            sep="\t", na="NA", row.names=FALSE, col.names=TRUE, quote=FALSE)


###############################################################################
#17. Supplementary table (published format)
#Reproduces SUPP_Table_Yeast_sps_paralogs_evo_trajectory_transcriptome.tsv:
#14 columns, one row per qualifying orthogroup, KEGG pathway titles as a
#semicolon-separated list.
###############################################################################
df_og2path <- df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY
df_og2path <- df_og2path[!is.na(df_og2path$ko_pathway) &
                           df_og2path$ko_pathway != "NA" &
                           df_og2path$ko_pathway != "" &
                           !is.na(df_og2path$pathway_title) &
                           df_og2path$pathway_title != "NA" &
                           df_og2path$pathway_title != "", ]

lst_path_titles <- split(df_og2path$pathway_title, df_og2path$Orthogroup)
v_kegg <- vapply(df_paralogs_fate$Orthogroup, function(og) {
  ti <- lst_path_titles[[og]]
  if (is.null(ti) || !length(ti)) return("NA")
  paste(unique(ti), collapse=";")            #order preserved, duplicates dropped
}, character(1))

df_supp <- data.frame(
  Orthogroup                      = df_paralogs_fate$Orthogroup,
  nb_species                      = df_paralogs_fate$nb_species,
  nb_unique_set_sizes             = df_paralogs_fate$nb_unique_set_sizes,
  ESD_rho                         = df_paralogs_fate$ESD_rho,
  ESD_vote                        = df_paralogs_fate$ESD_vote,
  slope_R2_same_rho               = df_paralogs_fate$slope_R2_same_rho,
  slope_R2_same_pval              = df_paralogs_fate$slope_R2_same_pval,
  R2slope_vote                    = df_paralogs_fate$R2slope_vote,
  markers_agree                   = df_paralogs_fate$markers_agree,
  min_abs_rho                     = df_paralogs_fate$min_abs_rho,
  Transcriptomic_trajectory_class = as.character(df_paralogs_fate$Transcriptomic_trajectory_class),
  unresolved_reason               = df_paralogs_fate$unresolved_reason,
  avg_truncation_score            = df_paralogs_fate$avg_truncation_score,
  KEGG_pathways                   = v_kegg,
  stringsAsFactors                = FALSE
)
write.table(df_supp,
            file = paste0(dir_fate_root,
                          "SUPP_Table_Yeast_sps_paralogs_evo_trajectory_transcriptome.tsv"),
            sep="\t", na="NA", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat("\nWrote supplementary table:", nrow(df_supp), "orthogroups x", ncol(df_supp), "columns\n")


###############################################################################
#18. Sensitivity analysis on the set-size threshold
#
#RATIONALE (reported in Methods): thresholds available to the sliding markers
#equal (distinct non-zero set sizes - 1), because the highest set size is
#dropped. The primary statistic requires >= 3 valid thresholds, so an
#orthogroup needs >= 4 distinct non-zero set sizes to be scored by it. Below
#that both markers fall back to a different per-species statistic. The
#threshold is therefore an implementation boundary, not a tuning choice.
#
#Two independent checks are reported:
#  (a) calls as a function of the minimum distinct set sizes required
#  (b) whether the class distribution at the lowest resolution is
#      distinguishable from what two random-sign markers would produce
###############################################################################
sensitivity_by_threshold <- function(df, thresholds = 3:7) {
  do.call(rbind, lapply(thresholds, function(k) {
    sub <- df[df$nb_unique_set_sizes >= k, ]
    tb  <- table(factor(sub$Transcriptomic_trajectory_class,
                        levels=c(CLASS_DIV, CLASS_COH, CLASS_UNR)))
    D <- as.integer(tb[CLASS_DIV]); C <- as.integer(tb[CLASS_COH]); U <- as.integer(tb[CLASS_UNR])
    data.frame(min_distinct_set_sizes = k,
               n_orthogroups = nrow(sub),
               n_thresholds_available = k - 1L,
               statistic_used = ifelse(k >= 4, "primary (sliding)", "fallback for some"),
               Divergence = D, Cohesion = C, Unresolved = U,
               pct_cohesion_of_called = if ((C+D) > 0) round(100*C/(C+D), 1) else NA_real_,
               ratio_cohesion_to_divergence = if (D > 0) round(C/D, 2) else NA_real_,
               pct_unresolved = round(100*U/nrow(sub), 1),
               stringsAsFactors = FALSE)
  }))
}
df_sensitivity <- sensitivity_by_threshold(df_paralogs_fate)
cat("\n=== SENSITIVITY: calls by minimum distinct (non-zero) set sizes ===\n")
print(df_sensitivity, row.names=FALSE)
write.table(df_sensitivity,
            file = paste0(dir_fate_root, "SUPP_Table_setsize_threshold_sensitivity.tsv"),
            sep="\t", na="NA", row.names=FALSE, col.names=TRUE, quote=FALSE)

#(b) random-sign null at the lowest resolution.
#Two markers with independent random signs agree 50% of the time: 25% both
#Divergence, 25% both Cohesion, 50% disagree. Test each stratum against that.
random_sign_test <- function(df, k) {
  sub <- df[df$nb_unique_set_sizes == k, ]
  if (!nrow(sub)) return(NULL)
  obs <- c(sum(sub$Transcriptomic_trajectory_class == CLASS_DIV),
           sum(sub$Transcriptomic_trajectory_class == CLASS_COH),
           sum(sub$Transcriptomic_trajectory_class == CLASS_UNR))
  ct <- suppressWarnings(chisq.test(obs, p = c(0.25, 0.25, 0.50)))
  data.frame(distinct_set_sizes = k, n = nrow(sub),
             pct_Divergence = round(100*obs[1]/sum(obs),1),
             pct_Cohesion   = round(100*obs[2]/sum(obs),1),
             pct_Unresolved = round(100*obs[3]/sum(obs),1),
             chisq = round(unname(ct$statistic),3),
             pval  = signif(ct$p.value, 3),
             stringsAsFactors = FALSE)
}
df_randomsign <- do.call(rbind, lapply(sort(unique(df_paralogs_fate$nb_unique_set_sizes)),
                                       function(k) random_sign_test(df_paralogs_fate, k)))
cat("\n=== Class distribution vs a random-sign null (25% / 25% / 50%) ===\n")
print(df_randomsign, row.names=FALSE)
write.table(df_randomsign,
            file = paste0(dir_fate_root, "SUPP_Table_random_sign_null_by_resolution.tsv"),
            sep="\t", na="NA", row.names=FALSE, col.names=TRUE, quote=FALSE)


###############################################################################
#19. FIG5: trajectory calls + resolution limit-
#Same two-panel layout as before; only the numbers change.
###############################################################################
COL_DIV  <- "#8E44AD"   #Divergence  (purple)
COL_COH  <- "#16A085"   #Cohesion    (teal)
COL_GREY <- "#7A7A7A"   #under-powered
COL_GOLD <- "#C9A227"   #strong discordant

n_tot <- nrow(df_paralogs_fate)
n_div <- sum(df_paralogs_fate$Transcriptomic_trajectory_class == CLASS_DIV)
n_coh <- sum(df_paralogs_fate$Transcriptomic_trajectory_class == CLASS_COH)
n_unr <- sum(df_paralogs_fate$Transcriptomic_trajectory_class == CLASS_UNR)
n_strong <- sum(df_paralogs_fate$unresolved_reason == "strong_discordant_signal", na.rm=TRUE)
n_weak   <- n_unr - n_strong
pct_strong <- if (n_unr > 0) round(100*n_strong/n_unr) else NA_integer_
pct_weak   <- if (n_unr > 0) 100 - pct_strong else NA_integer_

df_panelA <- rbind(
  data.frame(class=CLASS_DIV, part="called",           n=n_div,    stringsAsFactors=FALSE),
  data.frame(class=CLASS_COH, part="called",           n=n_coh,    stringsAsFactors=FALSE),
  data.frame(class=CLASS_UNR, part="under-powered",    n=n_weak,   stringsAsFactors=FALSE),
  data.frame(class=CLASS_UNR, part="strong discordant",n=n_strong, stringsAsFactors=FALSE)
)
df_panelA$prop <- df_panelA$n / n_tot
df_panelA$class <- factor(df_panelA$class, levels=c(CLASS_DIV, CLASS_COH, CLASS_UNR))
df_panelA$part  <- factor(df_panelA$part,
                          levels=c("strong discordant","under-powered","called"))
df_panelA$fill_col <- ifelse(df_panelA$part == "strong discordant", COL_GOLD,
                             ifelse(df_panelA$part == "under-powered",     COL_GREY,
                                    ifelse(df_panelA$class == CLASS_DIV,          COL_DIV, COL_COH)))

lab_div <- paste0("Transcriptomic\ndivergence\n(n=", n_div, ")")
lab_coh <- paste0("Transcriptomic\ncohesion\n(n=", n_coh, ")")
lab_unr <- paste0("Unresolved\n(markers are split)\n(n=", n_unr, ")")

pA <- ggplot(df_panelA, aes(x=class, y=prop, fill=fill_col, group=part)) +
  geom_col(width=0.72, colour=NA) +
  scale_fill_identity() +
  geom_text(aes(label=n), position=position_stack(vjust=0.5),
            colour="white", fontface="bold", size=5) +
  scale_x_discrete(labels=c(lab_div, lab_coh, lab_unr)) +
  scale_y_continuous(expand=expansion(mult=c(0, 0.08))) +
  labs(x="Transcriptomic trajectory followed by expanding orthogroups",
       y="Proportion of expanding orthogroups") +
  theme_bw(base_size=13) +
  theme(panel.grid.major.x=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(colour="grey30"),
        legend.position="none")

#Panel B: unresolved rate vs copy-number resolution
df_paralogs_fate$ss_bin <- cut(df_paralogs_fate$nb_unique_set_sizes,
                               breaks=c(2,3,4,5,6,Inf),
                               labels=c("3","4","5","6",">=7"), right=TRUE)
df_panelB <- do.call(rbind, lapply(levels(df_paralogs_fate$ss_bin), function(b) {
  sub <- df_paralogs_fate[!is.na(df_paralogs_fate$ss_bin) & df_paralogs_fate$ss_bin == b, ]
  data.frame(bin=b, n=nrow(sub),
             n_unres=sum(sub$Transcriptomic_trajectory_class == CLASS_UNR),
             stringsAsFactors=FALSE)
}))
df_panelB$pct_unres <- 100 * df_panelB$n_unres / df_panelB$n
df_panelB <- df_panelB[df_panelB$n > 0, ]          #prop.trend.test requires n > 0
df_panelB$bin <- factor(df_panelB$bin, levels=df_panelB$bin)

#Cochran-Armitage trend test on ordered bins (prop.trend.test = CA trend test).
#Two-sided by construction; halved for the directional (decreasing) hypothesis.
ca <- prop.trend.test(df_panelB$n_unres, df_panelB$n, score=seq_len(nrow(df_panelB)))
slope_dir <- coef(lm(df_panelB$pct_unres ~ seq_len(nrow(df_panelB))))[2]  #sign only
p_ca_one  <- if (slope_dir < 0) ca$p.value/2 else 1 - ca$p.value/2
#chi-square test of independence across bins (2 x k)
chi <- suppressWarnings(chisq.test(rbind(df_panelB$n_unres,
                                         df_panelB$n - df_panelB$n_unres)))
fmt_p <- function(p) formatC(p, format="e", digits=2)
lab_tests <- paste0("Cochran-Armitage decreasing trend test p = ", fmt_p(p_ca_one),
                    "\nChi-square test p = ", fmt_p(chi$p.value))

#endpoints of the least-squares trend across the ordered bins (visual guide only;
#the reported test is Cochran-Armitage, which uses the ordered ranks)
.xnum   <- seq_len(nrow(df_panelB))
.fit    <- lm(df_panelB$pct_unres ~ .xnum)
df_trend <- data.frame(x1 = 1, x2 = nrow(df_panelB),
                       y1 = unname(predict(.fit, data.frame(.xnum = 1))),
                       y2 = unname(predict(.fit, data.frame(.xnum = nrow(df_panelB)))))

pB <- ggplot(df_panelB, aes(x=bin, y=pct_unres)) +
  geom_col(fill=COL_GREY, width=0.72) +
  geom_text(aes(label=paste0(round(pct_unres), "%")), vjust=-0.6,
            fontface="bold", size=4.6) +
  geom_text(aes(y=3, label=paste0("n=", n)), colour="white",
            fontface="bold", size=4.2) +
  geom_segment(data=df_trend, aes(x=x1, xend=x2, y=y1, yend=y2),
               colour="#C0392B", linetype="dashed", linewidth=1.1,
               inherit.aes=FALSE) +
  annotate("label", x=3.6, y=max(df_panelB$pct_unres)*1.18, label=lab_tests,
           colour="#C0392B", fontface="bold", size=4.1, hjust=0.5,
           label.size=0.7, fill="white") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.22))) +
  labs(x="Distinct set sizes across species", y="% of orthogroups unresolved") +
  theme_bw(base_size=13) +
  theme(panel.grid.major.x=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(colour="grey30"))

fig5 <- ggarrange(pA, pB, ncol=2, labels=c("A","B"),
                  font.label=list(size=22, face="bold"), widths=c(1, 1))
if (requireNamespace("svglite", quietly=TRUE)) {
  ggsave(paste0(dir_fate_root, "FIG5_Expanding_orthogroups_trajectory_calls.svg"),
         fig5, width=15, height=6.2, dpi=600)
} else {                       #fall back to the base grDevices SVG device
  svg(paste0(dir_fate_root, "FIG5_Expanding_orthogroups_trajectory_calls.svg"),
      width=15, height=6.2); print(fig5); dev.off()
}
ggsave(paste0(dir_fate_root, "FIG5_Expanding_orthogroups_trajectory_calls.png"),
       fig5, width=15, height=6.2, dpi=600)
cat("\nWrote FIG5 (svg + png)\n")
cat("  Panel A:", n_div, "Divergence /", n_coh, "Cohesion /", n_unr, "Unresolved",
    "(", pct_weak, "% under-powered,", pct_strong, "% strong discordant )\n")
cat("  Panel B: Cochran-Armitage one-sided p =", fmt_p(p_ca_one),
    "; chi-square p =", fmt_p(chi$p.value), "\n")


###############################################################################
#20. Non-circular benchmark of candidate third markers
#
#The trajectory class is a DERIVED construct, so there is no external ground
#truth. Candidates are therefore scored on the orthogroups where the two
#class-defining markers AGREE: a target defined by those markers alone. The
#majority class on that set fixes the accuracy a non-informative constant would
#reach; a candidate must exceed it to carry independent information.
###############################################################################
df_agree <- df_paralogs_fate[df_paralogs_fate$markers_agree, ]
if (!nrow(df_agree)) stop("No orthogroups where both markers agree; cannot benchmark.")
target   <- as.character(df_agree$Transcriptomic_trajectory_class)
baseline <- max(mean(target == CLASS_COH), mean(target == CLASS_DIV))

#each candidate maps a row to a predicted class (or NA when not computable)
lst_candidates <- list(
  `R2 ratio (same vs different pathway)` =
    function(d) ifelse(is.na(d$slope_R2_ratio_rho), NA_character_,
                       ifelse(d$slope_R2_ratio_rho > 0, CLASS_COH, CLASS_DIV)),
  `Copy-number attrition` =
    function(d) ifelse(is.na(d$log_decay), NA_character_,
                       ifelse(d$log_decay, CLASS_DIV, CLASS_COH)),
  `Coefficient-of-variation slope` =
    function(d) ifelse(is.na(d$norm_CV_slope), NA_character_,
                       ifelse(d$norm_CV_slope < 0, CLASS_COH, CLASS_DIV)),
  `Partial dependence of expression on copy number` =
    function(d) ifelse(is.na(d$beta_setsize), NA_character_,
                       ifelse(d$beta_setsize > 0, CLASS_COH, CLASS_DIV)),
  `Nb species pairs substantially shortened` =
    function(d) ifelse(is.na(d$avg_nb_species_pairs_sign_trunc), NA_character_,
                       ifelse(d$avg_nb_species_pairs_sign_trunc > 0, CLASS_DIV, CLASS_COH)),
  `Gene-truncation score` =
    function(d) ifelse(is.na(d$avg_truncation_score), NA_character_,
                       ifelse(d$avg_truncation_score > 0.1, CLASS_DIV, CLASS_COH))
)

boot_lift <- function(pred, truth, n_boot = 2000) {
  ok <- !is.na(pred)
  if (sum(ok) < 10) return(c(NA_real_, NA_real_))
  p <- pred[ok]; t <- truth[ok]
  bl <- max(mean(t == CLASS_COH), mean(t == CLASS_DIV))
  lifts <- replicate(n_boot, {
    i <- sample(seq_along(p), replace=TRUE)
    mean(p[i] == t[i]) - max(mean(t[i] == CLASS_COH), mean(t[i] == CLASS_DIV))
  })
  unname(quantile(lifts, c(0.025, 0.975), na.rm=TRUE))
}

df_bench <- do.call(rbind, lapply(names(lst_candidates), function(nm) {
  pred <- lst_candidates[[nm]](df_agree)
  ok   <- !is.na(pred)
  acc  <- if (sum(ok)) mean(pred[ok] == target[ok]) else NA_real_
  ci   <- boot_lift(pred, target)
  #rank correlation with each primary marker, to flag redundancy
  num <- switch(nm,
                `R2 ratio (same vs different pathway)`            = df_agree$slope_R2_ratio_rho,
                `Copy-number attrition`                           = as.numeric(df_agree$log_decay),
                `Coefficient-of-variation slope`                  = df_agree$norm_CV_slope,
                `Partial dependence of expression on copy number` = df_agree$beta_setsize,
                `Nb species pairs substantially shortened`        = df_agree$avg_nb_species_pairs_sign_trunc,
                `Gene-truncation score`                           = df_agree$avg_truncation_score)
  rho_esd <- suppressWarnings(cor(num, df_agree$ESD_rho,
                                  method="spearman", use="pairwise.complete.obs"))
  rho_r2  <- suppressWarnings(cor(num, df_agree$slope_R2_same_rho,
                                  method="spearman", use="pairwise.complete.obs"))
  data.frame(candidate_marker = nm,
             n_evaluable      = sum(ok),
             agreement_pct    = round(100*acc, 1),
             constant_baseline_pct = round(100*baseline, 1),
             lift_pct         = round(100*(acc - baseline), 1),
             lift_CI95_low    = round(100*ci[1], 1),
             lift_CI95_high   = round(100*ci[2], 1),
             beats_constant   = !is.na(acc) && acc > baseline,
             rho_with_ESD     = round(rho_esd, 3),
             rho_with_R2slope = round(rho_r2, 3),
             stringsAsFactors = FALSE)
}))
df_bench <- df_bench[order(-df_bench$agreement_pct), ]
cat("\n=== Non-circular benchmark of candidate third markers ===\n")
cat("Target: the", nrow(df_agree), "orthogroups where both primary markers agree.\n")
cat("Constant baseline:", round(100*baseline,1), "%\n\n")
print(df_bench, row.names=FALSE)
write.table(df_bench,
            file = paste0(dir_fate_root, "SUPP_Table_candidate_marker_benchmark.tsv"),
            sep="\t", na="NA", row.names=FALSE, col.names=TRUE, quote=FALSE)


###############################################################################
#21. Numbers for the Methods / Results text
#Everything a downstream document needs, printed once and written to file.
###############################################################################
sink(paste0(dir_fate_root, "trajectory_analysis_key_numbers.txt"))
cat("TRAJECTORY ANALYSIS - KEY NUMBERS\n")
cat("=================================\n\n")
cat("Qualifying orthogroups:", n_tot, "\n")
cat("  Divergence :", n_div, "\n")
cat("  Cohesion   :", n_coh, "\n")
cat("  Unresolved :", n_unr, sprintf("(%.1f%%)", 100*n_unr/n_tot), "\n")
cat("  Called     :", n_div + n_coh,
    sprintf("-> %.1f%% cohesion, ratio %.2f:1\n", 100*n_coh/(n_div+n_coh), n_coh/max(n_div,1)))
cat("\nUnresolved partition:\n")
cat("  under-powered     :", n_weak,   sprintf("(%d%%)\n", pct_weak))
cat("  strong discordant :", n_strong, sprintf("(%d%%)\n", pct_strong))
cat("\nScoring path:\n")
print(table(df_paralogs_fate$scoring_path, useNA="ifany"))
cat("\nClass by scoring path:\n")
print(table(df_paralogs_fate$scoring_path,
            df_paralogs_fate$Transcriptomic_trajectory_class, useNA="ifany"))
cat("\nSet-size threshold sensitivity:\n")
print(df_sensitivity, row.names=FALSE)
cat("\nRandom-sign null by resolution:\n")
print(df_randomsign, row.names=FALSE)
cat("\nPanel B trend tests:\n")
cat("  Cochran-Armitage one-sided p =", fmt_p(p_ca_one), "\n")
cat("  Chi-square (2 x k) p         =", fmt_p(chi$p.value), "\n")
cat("\nCandidate-marker benchmark (baseline",
    sprintf("%.1f%%", 100*baseline), "):\n")
print(df_bench, row.names=FALSE)
cat("\nTruncation score by trajectory class (Kruskal-Wallis):\n")
df_kw <- df_paralogs_fate[!is.na(df_paralogs_fate$avg_truncation_score), ]
df_kw$Transcriptomic_trajectory_class <- droplevels(df_kw$Transcriptomic_trajectory_class)
if (nlevels(df_kw$Transcriptomic_trajectory_class) >= 2) {
  print(kruskal.test(avg_truncation_score ~ Transcriptomic_trajectory_class, data=df_kw))
} else cat("  <2 non-empty classes; test not run\n")
cat("\nCore vs accessory (nb_species == 16 vs < 16) among CALLED orthogroups:\n")
df_called2 <- df_paralogs_fate[df_paralogs_fate$Transcriptomic_trajectory_class != CLASS_UNR, ]
df_called2$compartment <- ifelse(df_called2$nb_species == length(v_sps_YPD),
                                 "core", "accessory")
tb_comp <- table(df_called2$compartment,
                 droplevels(df_called2$Transcriptomic_trajectory_class))
print(tb_comp)
if (all(dim(tb_comp) == c(2,2))) print(fisher.test(tb_comp)) else
  cat("  contingency table is not 2x2; Fisher test not run\n")
sink()
cat("\nWrote key numbers to trajectory_analysis_key_numbers.txt\n")
cat("\nAll outputs written to:", dir_fate_root, "\n")

###############################################################################
#Pie chart: orthogroup counts
#  - total orthogroups
#  - orthogroups WITH paralogs (>=2 copies in at least one species)
#  - orthogroups with paralogs whose SIZE INCREASES WITH TIME
#
#Because these three counts are NESTED (each is a subset of the previous), a
#pie of the raw counts would double-count. Instead the pie shows three
#MUTUALLY-EXCLUSIVE partitions that sum to the total:
#    (a) single-copy (no paralogs)
#    (b) has paralogs, NOT increasing with time
#    (c) has paralogs, increasing with time
#All three requested numbers are still reported (subtitle + console).
#
#Requires from the main script:
#  mtx_species_Orthogroups_set_size_YPD   rows=species (smpl_lbl), cols=OG; NA=absent
#  tip_distances                          named numeric (names = Genus_species), distance to root
#  output_workspace                       character path ending with "/"
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
})

#Parameters
MIN_SP_FOR_TREND <- 3        #min species present to test a size-vs-time trend
TREND_REQUIRE_SIG <- TRUE    #TRUE: Spearman rho>0 AND p<0.05; FALSE: rho>0 only
TREND_ALPHA <- 0.05

#1. Inputs
ss <- mtx_species_Orthogroups_set_size_YPD          #species x orthogroups
n_total <- ncol(ss)

#Per-orthogroup maximum copy number across species (NA-safe)
max_ss <- apply(ss, 2, function(x) { v <- x[!is.na(x)]; if (length(v)) max(v) else 0 })
has_paralogs <- max_ss > 1                          #>=2 copies in >=1 species

#2. Distance-to-root per species, aligned to matrix rows
smpl_lbl_to_tip_key <- function(sp) gsub(" ", "_", sp, fixed = TRUE)
root_dist <- tip_distances[smpl_lbl_to_tip_key(rownames(ss))]
names(root_dist) <- rownames(ss)

#3. Does an orthogroup increase in size with time?
#"time" = phylogenetic distance to root (more derived species = more time elapsed).
#Positive Spearman correlation between per-species set size and root distance.
increases_with_time <- function(og) {
  x <- ss[, og]
  keep <- !is.na(x) & x >= 1                        #species where OG is present
  if (sum(keep) < MIN_SP_FOR_TREND) return(FALSE)
  s <- x[keep]; d <- root_dist[keep]
  if (length(unique(s)) < 2 || length(unique(d)) < 2) return(FALSE)
  ct <- suppressWarnings(cor.test(d, s, method = "spearman", exact = FALSE))
  if (is.na(ct$estimate)) return(FALSE)
  if (TREND_REQUIRE_SIG) ct$estimate > 0 && ct$p.value < TREND_ALPHA
  else                   ct$estimate > 0
}

og_para  <- colnames(ss)[has_paralogs]
incr_flag <- vapply(og_para, increases_with_time, logical(1))

#4. Counts
n_with_paralogs <- length(og_para)
n_incr_time     <- sum(incr_flag)

n_single        <- n_total - n_with_paralogs
n_para_no_trend <- n_with_paralogs - n_incr_time

cat("Total orthogroups ............................. ", n_total, "\n")
cat("  with paralogs (>=2 copies in >=1 species) ... ", n_with_paralogs,
    sprintf(" (%.1f%%)\n", 100 * n_with_paralogs / n_total))
cat("    of which: size increases with time ........ ", n_incr_time,
    sprintf(" (%.1f%% of paralog OGs)\n", 100 * n_incr_time / max(n_with_paralogs, 1)))

#5. Pie data (mutually exclusive partitions that sum to n_total)
df_pie <- data.frame(
  category = factor(
    c("Single-copy (no paralogs)",
      "Paralogs, not increasing with time",
      "Paralogs, increasing with time"),
    levels = c("Single-copy (no paralogs)",
               "Paralogs, not increasing with time",
               "Paralogs, increasing with time")),
  count = c(n_single, n_para_no_trend, n_incr_time)
)
df_pie$pct   <- 100 * df_pie$count / sum(df_pie$count)
df_pie$label <- sprintf("%s\n%d (%.1f%%)", df_pie$category, df_pie$count, df_pie$pct)

#colorblind-friendly (Paul Tol)
pal <- c("Single-copy (no paralogs)"          = "#BBBBBB",
         "Paralogs, not increasing with time" = "#4477AA",
         "Paralogs, increasing with time"     = "#CC6677")

#6. Pie chart
sub <- sprintf("Total = %d orthogroups  |  with paralogs = %d  |  paralogs increasing with time = %d",
               n_total, n_with_paralogs, n_incr_time)

p_pie <- ggplot(df_pie, aes(x = "", y = count, fill = category)) +
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = pal, name = NULL) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", count, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.6, color = "white", fontface = "bold") +
  labs(title = "Orthogroups, paralog content, and copy-number expansion over time",
       subtitle = sub) +
  theme_void(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5,
                                     margin = margin(b = 6)),
        legend.position = "right",
        legend.text   = element_text(size = 9))

#7. Save (600 dpi PNG + SVG)
ggsave("Piechart_orthogroup_paralog_counts.png", p_pie, path = output_workspace,
       width = 22, height = 16, units = "cm", dpi = 600)
ggsave("Piechart_orthogroup_paralog_counts.svg", p_pie, path = output_workspace,
       width = 22 / 2.54, height = 16 / 2.54, device = svg)

cat("\nSaved: Piechart_orthogroup_paralog_counts.{png,svg} to", output_workspace, "\n")

################Show that pairs of core and accessory orthogroups have more correlation within a pathway than across different pathways
core_ogs <- intersect(v_lst_core_Orthogroups_YPD,  colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD))
acc_ogs  <- intersect(v_lst_accessory_Orthogroups_YPD, colnames(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD))
MIN_BOTH_EXPRESSED <- 3
media <- "YPD"

df_og_path <- unique(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY[
  !is.na(df_yeasts_pantranscriptome_WITH_DUPLICATES_FROM_KEGG_PATHWAY$ko_pathway),
  c("Orthogroup","ko_pathway")])
lst_og_pathways <- split(as.character(df_og_path$ko_pathway), df_og_path$Orthogroup)

core_ogs <- core_ogs[core_ogs %in% names(lst_og_pathways)]
acc_ogs  <- acc_ogs [acc_ogs  %in% names(lst_og_pathways)]
all_ogs <- c(core_ogs, acc_ogs)                    #disjoint (core vs accessory)

#same-pathway indicator via a binary membership matrix (fast)
all_paths <- sort(unique(unlist(lst_og_pathways[all_ogs])))
M <- matrix(0L, nrow = length(all_ogs), ncol = length(all_paths),
            dimnames = list(all_ogs, all_paths))
for (i in seq_along(all_ogs)) {
  pw <- lst_og_pathways[[all_ogs[i]]]; if (length(pw)) M[i, pw] <- 1L
}
ci <- seq_along(core_ogs); ai <- length(core_ogs) + seq_along(acc_ogs)
share_mat  <- (M[ci, , drop=FALSE] %*% t(M[ai, , drop=FALSE])) > 0   #n_core x n_acc
within_vec <- as.vector(share_mat)

#R^2 over ALL species (fill=0 kept)
R2_all <- (cor(mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[, core_ogs, drop=FALSE],
               mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[, acc_ogs,  drop=FALSE]))^2

#R^2 over species where BOTH are expressed (0 -> NA)
Ec <- mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[, core_ogs, drop=FALSE]; Ec[Ec == 0] <- NA
Ea <- mtx_log2_1p_species_aggregate_gene_expr_FPKM_YPD[, acc_ogs,  drop=FALSE]; Ea[Ea == 0] <- NA
R_both <- suppressWarnings(cor(Ec, Ea, use = "pairwise.complete.obs"))
n_both <- crossprod(+!is.na(Ec), +!is.na(Ea))        #core x acc: #both-expressed species
R_both[n_both < MIN_BOTH_EXPRESSED] <- NA            #need >=3 both-expressed species
R2_both <- R_both^2

#assemble long data frame
lab_all  <- "All species (0 kept)"
lab_both <- sprintf("Both expressed (>=%d sp.)", MIN_BOTH_EXPRESSED)
df_pairs <- rbind(
  data.frame(R2 = as.vector(R2_all),  within = within_vec, metric = lab_all,  stringsAsFactors = FALSE),
  data.frame(R2 = as.vector(R2_both), within = within_vec, metric = lab_both, stringsAsFactors = FALSE))
df_pairs$group  <- factor(ifelse(df_pairs$within, "Within pathway", "Between pathways"),
                          levels = c("Within pathway", "Between pathways"))
df_pairs$metric <- factor(df_pairs$metric, levels = c(lab_all, lab_both))
df_pairs <- df_pairs[!is.na(df_pairs$R2), ]

#per-metric Wilcoxon + label-permutation p
r2a <- as.vector(R2_all); r2b <- as.vector(R2_both)
report_metric <- function(vec, name) {
  wi <- vec[within_vec]; be <- vec[!within_vec]
  wi <- wi[!is.na(wi)]; be <- be[!is.na(be)]
  pw <- suppressWarnings(wilcox.test(wi, be, alternative = "greater")$p.value)
  message(sprintf("[%s] %-26s within n=%d med=%.4f | between n=%d med=%.4f | Wilcoxon p=%.3g",
                  media, name, length(wi), median(wi), length(be), median(be), pw))
  invisible(median(wi) - median(be))
}
obs_a <- report_metric(r2a, lab_all)
obs_b <- report_metric(r2b, lab_both)

set.seed(1)
pa <- pb <- numeric(N_PERM)
for (p in seq_len(N_PERM)) {
  Mp <- M[sample(nrow(M)), , drop=FALSE]                      #shuffle pathway labels across OGs
  Sp <- as.vector((Mp[ci, , drop=FALSE] %*% t(Mp[ai, , drop=FALSE])) > 0)
  pa[p] <- median(r2a[Sp], na.rm=TRUE) - median(r2a[!Sp], na.rm=TRUE)
  pb[p] <- median(r2b[Sp], na.rm=TRUE) - median(r2b[!Sp], na.rm=TRUE)
}
fmtp <- function(pv) if (pv <= 1/(N_PERM+1)) paste0("<", signif(1/(N_PERM+1),2)) else signif(pv,3)
p_perm_a <- (1 + sum(pa >= obs_a, na.rm=TRUE)) / (N_PERM + 1)
p_perm_b <- (1 + sum(pb >= obs_b, na.rm=TRUE)) / (N_PERM + 1)
message(sprintf("[%s] permutation p (within>between): all-species %s | both-expressed %s",
                media, fmtp(p_perm_a), fmtp(p_perm_b)))

#faceted violin+boxplot with stat_compare_means()
comparisons <- list(c("Within pathway", "Between pathways"))
gg <- ggplot(df_pairs, aes(x = group, y = R2)) +
  geom_violin(fill = "grey") +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA) +
  stat_compare_means(comparisons = comparisons, method = "wilcox",
                     method.args = list(alternative = "greater")) +
  facet_wrap(~ metric, scales = "free_y") +
  xlab("") + ylab(expression("Core-accessory pairwise expression R"^2)) +
  ggtitle(sprintf("Core-accessory co-expression within vs between KEGG pathways (%s)", media)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 11),
        strip.text = element_text(size = 12), plot.title = element_text(size = 12))

ggsave(paste0("Boxplot_CoreAccessory_within_vs_between_pathway_R2_", media, ".png"),
       gg, path = output_workspace, width = 20, height = 15, units = "cm", dpi = 600)
ggsave(paste0("Boxplot_CoreAccessory_within_vs_between_pathway_R2_", media, ".svg"),
       gg, path = output_workspace, width = 20/2.54, height = 15/2.54, device = svg)
################
#Find the number of expanding orthogroups
anc <- df_all_pairs$sp_ancient
der <- df_all_pairs$sp_derived

A <- mtx_species_Orthogroups_set_size_YPD[anc, , drop = FALSE]   #pairs x orthogroups, ancestral member
D <- mtx_species_Orthogroups_set_size_YPD[der, , drop = FALSE]   #pairs x orthogroups, derived member

#same conditions as log_decay_test: both species must carry >=1 copy, and ss_d - ss_a > 0
gain <- (!is.na(A) & A >= 1) & (!is.na(D) & D >= 1) & ((D - A) > 0)

nb_sp_pairs_in_which_each_orthogroup_has_expanded <- n_expanding_pairs <- colSums(gain)                    #per orthogroup
nb_expanding_orthogroups <- sum(n_expanding_pairs > 0) 

################Within vs between pathway R2###################
#add mean_R2 stats
v_ogs_R2 <- intersect(v_all_orthogroups_with_pathway, colnames(mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD))
cat("Computing mean R2 for", length(v_ogs_R2), "orthogroups with pathway annotation",
    "(df_paralogs_fate had", nrow(df_paralogs_fate), ")\n")

pick <- function(x, a, b) {
  if (!is.list(x)) return(NA_real_)
  v <- if (!is.null(x[[a]])) x[[a]] else x[[b]]
  if (is.null(v) || !length(v) || !is.finite(v)) NA_real_ else as.numeric(v)
}

df_R2_means_all <- do.call(rbind, lapply(seq_along(v_ogs_R2), function(i) {
  og <- v_ogs_R2[i]
  if (i %% 250 == 0) cat("  ", i, "/", length(v_ogs_R2), "\n")
  r2 <- tryCatch(compute_og_R2_with_pathway(og, mtx_species_Orthogroups_log2_1p_mean_standardized_FPKM_YPD), error = function(e) NULL)
  data.frame(Orthogroup   = og,
             mean_R2_same = pick(r2, "mean_R2_same", "avg_R2_same"),
             mean_R2_diff = pick(r2, "mean_R2_diff", "avg_R2_diff"),
             stringsAsFactors = FALSE)
}))

cat("Non-NA mean_R2_same:", sum(!is.na(df_R2_means_all$mean_R2_same)),
    "/ mean_R2_diff:",      sum(!is.na(df_R2_means_all$mean_R2_diff)),
    "of", nrow(df_R2_means_all), "\n")

df_src <- df_R2_means_all
keep   <- !is.na(df_src$mean_R2_same) & !is.na(df_src$mean_R2_diff)
df_R2_paired <- data.frame(Orthogroup = df_src$Orthogroup[keep],
                           R2_same    = df_src$mean_R2_same[keep],
                           R2_diff    = df_src$mean_R2_diff[keep],
                           stringsAsFactors = FALSE)
df_R2_paired <- df_R2_paired[order(df_R2_paired$Orthogroup), ]
n_og <- nrow(df_R2_paired)
cat("Orthogroups in within-vs-between R2 comparison:", n_og, "\n")

df_R2_long <- rbind(
  data.frame(Orthogroup = df_R2_paired$Orthogroup, mean_R2 = df_R2_paired$R2_same,
             pathway_relation = rep("Same pathway", n_og),      stringsAsFactors = FALSE),
  data.frame(Orthogroup = df_R2_paired$Orthogroup, mean_R2 = df_R2_paired$R2_diff,
             pathway_relation = rep("Different pathway", n_og), stringsAsFactors = FALSE))
df_R2_long$pathway_relation <- factor(df_R2_long$pathway_relation,
                                      levels = c("Same pathway", "Different pathway"))

#boxplot mean R2 same pathway vs different pathways
ggboxplot(df_R2_long, x = "pathway_relation", y = "mean_R2") +
  stat_compare_means(
    method = "wilcox.test",
    method.args = list(paired = TRUE, alternative = "greater"), label = "p.format",
    label.y = 1  #Adjust position of the p-value label
  ) + xlab("Are the orthogroups in the same pathway?") +
  ylab(expression("Orthogroup mean expression " * R^2)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme(legend.position = "none")
ggsave(filename = "Orthogroup_mean_expr_Rsq_within_vs_between_pathway_in_YPD.png", path=output_workspace, width = 15/2.54, height = 10/2.54)
ggsave(filename = "Orthogroup_mean_expr_Rsq_within_vs_between_pathway_in_YPD.svg", path=output_workspace, width = 15/2.54, height = 10/2.54,device=svg)


######################################################################


#save session
session::save.session(file = paste0(output_workspace,"Yeast_Pantranscriptomics_session_",gsub(pattern = ":",replacement = "_",x = gsub(pattern = " ",replacement = "_",x = date())),"_RSession.Rda"))
