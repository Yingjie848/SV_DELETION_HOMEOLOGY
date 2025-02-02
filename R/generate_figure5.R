# Generate plots/data for the figure5 and supplementary figures in the manuscript
# Here are some notes on the parameters used in the analysis:
#  1) we use 80% similarity, 30bp alignment length
#  2) we use the most similar alignment
#  3) we use the deletions with length 1kb-100kb in the main results, and deletions with length 1kb-10kb, 10kb-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb from Serena's paper in the supplementary results
#  4) we use Serena's data as main figure and Marcin's BOPP data as supplementary figure
#  5) we compare BRCA2 and control in the main results, and BRCA1, BRCA2, control in the supplementary results
# Note: propDelHome: proportion of deletions with homeology

rm(list=ls())

source("lib/lib_generate_data_for_figure5.R")

# set directory for the figures
figDir <- "figures"; dir.create(figDir, showWarnings = FALSE, recursive = TRUE)


## -----------------------------------------------------------------------------------------------
# Serena's 560 breast tumors

# set the directory of the pipeline output
dataDir="output/serena_SV_deletions_ssa_blast_lr100bp"; dir.create(dataDir)

# Serena01. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-100kb
make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena01_prop_of_dels_with_homeology_in_brca2_vs_control"))

# Serena02. Histogram of deletion lengths: separate histograms for control and BRCA2 (keeping the Y-axis the same for both control and BRCA2) 
make_histogram_of_deletion_lengths_for_brca2_and_control(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"Serena02_histogram_of_deletion_lengths_for_brca2_and_control/updated"))

# Serena03. Histogram of homeology lengths in BRCA2
make_histogram_of_homeology_length_in_brca2(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"Serena03_histogram_of_homeology_lengths_in_brca2/updated"))

# Serena04. Location of breakpoints (Serena): pie charts for breakpoints of deletions with homology and all deletions – 2 charts each for indels and SVs, and pie chart of frequency of repeat regions in the genome
make_breakpoint_located_on_repeats_plots(
    dataPathAllDeletions=file.path(dataDir,"output_deletions_1KB-100KB_final/results_repeats/all_deletions/repeats_summary.tsv"),
    dataPathDelHomeology=file.path(dataDir,"output_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.tsv"),
    dataPathDelHomeologyRepclass=file.path(dataDir,"output_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.repClass.tsv"),
    outDir=file.path(figDir,"Serena04_breakpoint_located_on_repeats_plots"))

# Serena05. Boxplot of frequency of deletions with homeology in HRDetect high vs low (Serena, excluding BRCA1 and BRCA2) for indels and SVs
make_boxplot_of_proportion_of_deletions_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena05_prop_of_dels_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2"))

# Serena06. POLQ expression in HRDetect high vs low (Serena, BRCA1 and BRCA2 excluded)
make_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena06_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2"))

# Serena07. RAD52 expression in HRDetect high vs low (Serena, BRCA1 and BRCA2 excluded)
make_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena07_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2"))

# Serena08. POLQ expression in BRCA2 and control
make_POLQ_expression_in_BRCA2_and_control(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena08_POLQ_expression_in_BRCA2_and_control"))

# Serena09. RAD52 expression in BRCA2 and control
make_RAD52_expression_in_BRCA2_and_control(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena09_RAD52_expression_in_BRCA2_and_control"))

# Serena10. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-10kb, 10kb-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control_for_different_deletion_lengths(
    dataPath1KB_10KB=file.path(dataDir,"output_deletions_1KBto10KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath10KB_100KB=file.path(dataDir,"output_deletions_10KBto100KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath100KB_1MB=file.path(dataDir,"output_deletions_100KBto1MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath1MB_10MB=file.path(dataDir,"output_deletions_1MBto10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath10MB=file.path(dataDir,"output_deletions_above_10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena10_prop_of_dels_with_homeology_in_brca2_vs_control_for_different_deletion_lengths"))


# Serena11. make piechart for proportion of repeat regions in the genome
make_piechart_for_proportion_of_repeat_regions_in_genome(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/results_repeats/genome_repeats/repeats_summary_in_genome.tsv"),
    outDir=file.path(figDir,"Serena11_piechart_for_proportion_of_repeat_regions_in_genome"))

# Serena12. make scatterplot between deletion and homeology length
make_scatterplot_between_deletion_and_homeology_length(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"Serena12_scatterplot_between_deletion_and_homeology_length"))

# Serena13. make gene expression heatmap
#make_gene_expression_heatmap(
#    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
#    outDir=file.path(figDir,"Serena13_gene_expression_heatmap")
#)

#make_gene_expression_heatmap(
#    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
#    outDir=file.path(figDir,"Serena13_gene_expression_heatmap/trimmed_values")
#)

make_gene_expression_heatmap(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena13_gene_expression_heatmap/trimmed_values_20240628")
)



## -----------------------------------------------------------------------------------------------
# Marcin's BOPP data

# set the directory of the pipeline output
dataDir="output/marcin_allTumor_lr100_blast"; dir.create(dataDir)


# BOPP01. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-100kb
make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control(
    dataPath=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"BOPP01_prop_of_dels_with_homeology_in_brca2_vs_control"))


# BOPP02. Histogram of deletion lengths: separate histograms for control and BRCA2 (keeping the Y-axis the same for both control and BRCA2) 
make_histogram_of_deletion_lengths_for_brca2_and_control(
    dataPath=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"BOPP02_histogram_of_deletion_lengths_for_brca2_and_control/updated"))

# BOPP03. Histogram of homeology lengths in BRCA2
make_histogram_of_homeology_length_in_brca2(
    dataPath=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"BOPP03_histogram_of_homeology_lengths_in_brca2/updated"))

# BOPP04. Location of breakpoints: pie charts for breakpoints of deletions with homology and all deletions – 2 charts each for indels and SVs, and pie chart of frequency of repeat regions in the genome
make_breakpoint_located_on_repeats_plots(
    dataPathAllDeletions=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/results_repeats/all_deletions/repeats_summary.tsv"),
    dataPathDelHomeology=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.tsv"),
    dataPathDelHomeologyRepclass=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.repClass.tsv"),
    outDir=file.path(figDir,"BOPP04_breakpoint_located_on_repeats_plots"))

# BOPP05. Boxplot of frequency of deletions with homeology in HRDetect high vs low (excluding BRCA1 and BRCA2) for indels and SVs
make_boxplot_of_proportion_of_deletions_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2(
    dataPath=file.path(dataDir,"BOPP_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"BOPP05_prop_of_dels_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2"))

# BOPP06. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-10kb, 10kb-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control_for_different_deletion_lengths(
    dataPath1KB_10KB=file.path(dataDir,"BOPP_deletions_1KBto10KB/ssa_events_samples.tsv"),
    dataPath10KB_100KB=file.path(dataDir,"BOPP_deletions_10KBto100KB/ssa_events_samples.tsv"),
    dataPath100KB_1MB=file.path(dataDir,"BOPP_deletions_100KBto1MB/ssa_events_samples.tsv"),
    dataPath1MB_10MB=file.path(dataDir,"BOPP_deletions_1MBto10MB/ssa_events_samples.tsv"),
    dataPath10MB=file.path(dataDir,"BOPP_deletions_above_10MB/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"BOPP06_prop_of_dels_with_homeology_in_brca2_vs_control_for_different_deletion_lengths"))

# BOPP07. Boxplot of frequency of deletions with homeology in HRDetect high vs low (BOPP, excluding BRCA1 and BRCA2)
make_boxplots_for_BOPP_tumor_types(
    dataPath=file.path(dataDir,"/BOPP_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    hrd_tbl_path="input/setton_hadi_choo_2023/hrd-supp-table_ZC.rds",
    outDir=file.path(figDir,"BOPP07_boxplots_for_BOPP_tumor_types")
)


## -----------------------------------------------------------------------------------------------
# summary of the results

# Serena01. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-100kb
d <- fread("figures/Serena01_prop_of_dels_with_homeology_in_brca2_vs_control/boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control.tsv")
d %>% group_by(allelic_status_brca1_brca2) %>% summarise(n=n())
#allelic_status_brca1_brca2     n
#1 BRCA2                         30
#2 control                      335

# Serena05. Boxplot of frequency of deletions with homeology in HRDetect high vs low (Serena, excluding BRCA1 and BRCA2) for indels and SVs
d <- fread("figures/Serena05_prop_of_dels_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2/boxplot_for_proportion_of_deletions_with_homeology_in_HRDetect_high_low_excluding_BRCA1_BRCA2.tsv")
d %>% group_by(HRDetect_group) %>% summarise(n=n())
#HRDetect_group     n
#1 HRDetect-high     44
#2 HRDetect-low     330


data <- fread(file.path("output/serena_SV_deletions_ssa_blast_lr100bp/output_deletions_1KB-100KB_final/ssa_events_samples.tsv")) %>% dplyr::filter(deletions != 0)
table(data$allelic_status_brca1_brca2,data$HRDetect_group)
#          HRDectect-intermediate HRDetect-high HRDetect-low
#  BRCA1                        6            41            0
#  BRCA2                        3            25            2
#  control                     12             5          318
#  other                        4            39           12


# BOPP01. Boxplot for frequency of deletions + homeology in BRCA2 vs control from Serena dataset 1kb-100kb
d <- fread("figures/BOPP01_prop_of_dels_with_homeology_in_brca2_vs_control/boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control.tsv")
d %>% group_by(allelic_status_brca1_brca2) %>% summarise(n=n())
#allelic_status_brca1_brca2     n
#1 BRCA2                        127
#2 control                      680

# BOPP05. Boxplot of frequency of deletions with homeology in HRDetect high vs low (excluding BRCA1 and BRCA2) for indels and SVs
d <- fread("figures/BOPP05_prop_of_dels_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2/boxplot_for_proportion_of_deletions_with_homeology_in_HRDetect_high_low_excluding_BRCA1_BRCA2.tsv")
d %>% group_by(HRDetect_group) %>% summarise(n=n())
#HRDetect_group     n
#1 HRDetect-high     16
#2 HRDetect-low     652

data <- fread(file.path("output/marcin_allTumor_lr100_blast/output_deletions_1KB-100KB_final/ssa_events_samples.tsv")) %>% dplyr::filter(deletions != 0)
table(data$allelic_status_brca1_brca2,data$HRDetect_group)
#          HRDectect-intermediate HRDetect-high HRDetect-low
#  BRCA1                        5            68           20
#  BRCA2                        6           125           19
#  control                     81            31         2954

# BOPP07. Boxplot of frequency of deletions with homeology in HRDetect high vs low (BOPP, excluding BRCA1 and BRCA2)
d <- fread("figures/BOPP07_boxplots_for_BOPP_tumor_types/sample_propDelHome.tsv")
d %>% group_by(tumor_type_final,allelic_status_brca1_brca2) %>% summarise(n=n())
#  tumor_type_final allelic_status_brca1_brca2     n
#1 BRCA             BRCA2                         66
#2 BRCA             control                      187
#3 OV               BRCA2                         16
#4 OV               control                       21
#5 PACA             BRCA2                          9
#6 PACA             control                      159
#7 PRAD             BRCA2                         36
#8 PRAD             control                      313