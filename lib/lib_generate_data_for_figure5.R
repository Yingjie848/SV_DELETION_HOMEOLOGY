
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(ggpubr)
library(tibble)


make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena01_prop_of_dels_with_homeology_in_brca2_vs_control")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # make boxplot
    # set ymax = 0.5
    group = "allelic_status_brca1_brca2"
    comparisons = list(c('BRCA2','control'))
    
    dplot <- data %>% 
        dplyr::filter(deletions != 0)
    # save brca1, brca2, control
    dplot %>% fwrite(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_brca1_brca2_vs_control.tsv"), sep="\t")
    
    # save brca2, control
    dplot <- dplot %>%
        dplyr::filter(get(group) %in% unlist(comparisons))
    dplot %>% fwrite(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control.tsv"), sep="\t")

    p <- dplot %>%
        ggplot(aes(x = get(group), y = homeology_rate)) +
        geom_boxplot() +
        theme_prism() +
        xlab("") +
        ylab("Proportion of deletions with homeology") +
        stat_compare_means(
            comparisons = comparisons,
            method = "wilcox.test", method.args = list("exact" = FALSE),
            label.y=0.4
        ) +
        coord_cartesian(ylim = c(0, 0.5))
        # theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
    ggsave(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control.ymax0.5.pdf"), width = 4, height = 5)

}


make_histogram_of_deletion_lengths_for_brca2_and_control <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"Serena02_histogram_of_deletion_lengths_for_brca2_and_control")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # prepare data    
    group = "allelic_status_brca1_brca2"
    columns = c("SAMPLE.TUMOR","CHROM","start_position","end_position","del_length","alignment_length","allelic_status_brca1_brca2")

    # make histogram, normalized in each group
    dplot <- data %>% dplyr::select(any_of(columns)) %>%
        dplyr::filter(get(group) %in% c('BRCA2','control'))

    p <- dplot %>% ggplot(aes(x=del_length,y=stat(density*width))) +  # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
        geom_histogram(color="black", binwidth=0.2) + 
        theme_prism() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + 
        ggtitle("Deletion length for deletions with homeology") + scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len_histgram_in_BRCA2_vs_control.log10.pdf"),width=7,height=4)

    p <- dplot %>% ggplot(aes(x=del_length,y=after_stat(count/sum(count)))) +  # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
        geom_histogram(color="black", binwidth=0.2) + 
        theme_prism() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + 
        geom_rug(aes(y=NULL)) +
        ggtitle("Deletion length for deletions with homeology") + scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len_histgram_in_BRCA2_vs_control.log10_rug.pdf"),width=7,height=4)

    p <- dplot %>% ggplot(aes(x=del_length,y=after_stat(count/sum(count)))) +  # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
        geom_histogram(color="black", binwidth=5000) + 
        theme_prism() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + 
        geom_rug(aes(y=NULL)) +
        ggtitle("Deletion length for deletions with homeology") + #scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len_histgram_in_BRCA2_vs_control.updated1.pdf"),width=7,height=4)

    p <- dplot %>% ggplot(aes(x=del_length,y=after_stat(count/sum(count)))) +  # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
        geom_histogram(color="black", binwidth=100) + 
        theme_prism() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + 
        geom_rug(aes(y=NULL)) +
        ggtitle("Deletion length for deletions with homeology") + #scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len_histgram_in_BRCA2_vs_control.updated2.pdf"),width=7,height=4)

    dplot %>% fwrite(paste0(outDir, "/deletion_len_histgram_in_BRCA2_vs_control.tsv"), sep="\t")

}


make_histogram_of_homeology_length_in_brca2 <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena03_histogram_of_homeology_lengths_in_brca2")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # normalized in each group
    group = "allelic_status_brca1_brca2"
    columns = c("SAMPLE.TUMOR","CHROM","start_position","end_position","del_length","alignment_length","allelic_status_brca1_brca2")

    # make histogram, normalized in each group
    dplot <- data %>% dplyr::select(any_of(columns)) %>%
        dplyr::filter(get(group) %in% 'BRCA2')

    p <- dplot %>% ggplot(aes(x=alignment_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black",binwidth=10) + 
        theme_prism() + xlab("Homeology length (bp)") + ylab("Proportion of deletions") + 
        geom_rug(aes(y=NULL)) +
        ggtitle("Homeology length")
    ggsave(paste0(outDir,"/homeology_len_histgram_BRCA2.pdf"),width=4,height=4)

    dplot %>% fwrite(paste0(outDir, "/homeology_len_histgram_BRCA2.tsv"), sep="\t")


    ## plot control group
    dplot <- data %>% dplyr::select(any_of(columns)) %>%
        dplyr::filter(get(group) %in% 'control')

    p <- dplot %>% ggplot(aes(x=alignment_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black",binwidth=10) + 
        theme_prism() + xlab("Homeology length (bp)") + ylab("Proportion of deletions") + 
        geom_rug(aes(y=NULL)) +
        ggtitle("Homeology length")
    ggsave(paste0(outDir,"/homeology_len_histgram_control.pdf"),width=4,height=4)

    dplot %>% fwrite(paste0(outDir, "/homeology_len_histgram_control.tsv"), sep="\t")

}


make_scatterplot_between_deletion_and_homeology_length <- function(
    dataPath=file.path(dataDir,"output_deletions_1KB-100KB_final/ssa_events_candidates.tsv"),
    outDir=file.path(figDir,"Serena12_scatterplot_between_deletion_and_homeology_length")){

    
    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # prepare data
    group = "allelic_status_brca1_brca2"
    columns = c("SAMPLE.TUMOR","CHROM","start_position","end_position","del_length","alignment_length","allelic_status_brca1_brca2")

    # for BRCA2 and control tumors
    dplot <- data %>% dplyr::select(any_of(columns)) %>%
        dplyr::filter(get(group) %in% c('BRCA2','control'))

    p <- dplot %>% ggplot(aes(x=del_length,y=alignment_length)) + 
        geom_point() + theme_prism() + xlab("Deletion length (bp)") + ylab("Homeology length (bp)") + 
        stat_cor(color='red') +
        ggtitle("Homeology length vs Deletion length") + scale_x_log10() + scale_y_log10()

    ggsave(paste0(outDir,"/scatterplot_between_deletion_and_homeology_length.BRCA2_control.pdf"),p,width=5,height=5)

    dplot %>% fwrite(paste0(outDir, "/scatterplot_between_deletion_and_homeology_length.BRCA2_control.tsv"), sep="\t")

    # for BRCA1, BRCA2, control tumors
    dplot <- data %>% dplyr::select(any_of(columns)) %>%
        dplyr::filter(get(group) %in% c('BRCA1','BRCA2','control'))

    p <- dplot %>% ggplot(aes(x=del_length,y=alignment_length)) + 
        geom_point() + theme_prism() + xlab("Deletion length (bp)") + ylab("Homeology length (bp)") + 
        stat_cor(color='red') +
        ggtitle("Homeology length vs Deletion length") + scale_x_log10() + scale_y_log10()

    ggsave(paste0(outDir,"/scatterplot_between_deletion_and_homeology_length.BRCA1_BRCA2_control.pdf"),p,width=5,height=5)

    dplot %>% fwrite(paste0(outDir, "/scatterplot_between_deletion_and_homeology_length.BRCA1_BRCA2_control.tsv"), sep="\t")

}

make_breakpoint_located_on_repeats_plots <- function(
    dataPathAllDeletions=file.path(dataDir,"/output_deletions_1KB-100KB_final/results_repeats/all_deletions/repeats_summary.tsv"),
    dataPathDelHomeology=file.path(dataDir,"/output_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.tsv"),
    dataPathDelHomeologyRepclass=file.path(dataDir,"/output_deletions_1KB-100KB_final/results_repeats/deletions_with_homeology/repeats_summary.repClass.tsv"),
    outDir=file.path(figDir,"Serena04_breakpoint_located_on_repeats_plots")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    dataAllDeletions <- fread(dataPathAllDeletions)
    dataDelHomeology <- fread(dataPathDelHomeology)
    dataDelHomeologyRepclass <- fread(dataPathDelHomeologyRepclass)

    # plot

    # make pie chart for all deletions
    dplot <- dataAllDeletions %>%
        pivot_longer(cols=c('pct_Both','pct_One','pct_No'),names_to="repeats_location",values_to="pct_repeats") %>%
        dplyr::mutate(repeats_location=case_when(
            repeats_location=="pct_Both" ~ "Both breakpoints",
            repeats_location=="pct_One" ~ "One breakpoint",
            repeats_location=="pct_No" ~ "No repeats",
        )) %>%
        dplyr::mutate(repeats_location=factor(repeats_location,levels=c("Both breakpoints","One breakpoint","No repeats")))

    p <- dplot %>%
        ggplot(aes(x="",y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", width=1) + 
            coord_polar("y", start=0) + 
            geom_text(aes(label = paste0(round(pct_repeats,0),"%")), size = 6, color = "black",
                position = position_stack(vjust = 0.5)) +
            theme_void() + scale_fill_hue(name="Repeats Location") +
            ggtitle("Deletions around repeats (both/one/no breakpoints)")
    ggsave(paste0(outDir,"/all_deletions.piechart.pdf"),width=5,height=5)

    dplot %>% fwrite(paste0(outDir, "/all_deletions.piechart.tsv"), sep="\t")

    # make pie chart for deletions with homeology
    dplot <- dataDelHomeology %>%
        pivot_longer(cols=c('pct_Both','pct_One','pct_No'),names_to="repeats_location",values_to="pct_repeats") %>%
        dplyr::mutate(repeats_location=case_when(
            repeats_location=="pct_Both" ~ "Both breakpoints",
            repeats_location=="pct_One" ~ "One breakpoint",
            repeats_location=="pct_No" ~ "No repeats",
        )) %>%
        dplyr::mutate(repeats_location=factor(repeats_location,levels=c("Both breakpoints","One breakpoint","No repeats")))

    p <- dplot %>%
        ggplot(aes(x="",y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", width=1) + 
            coord_polar("y", start=0) + 
            geom_text(aes(label = paste0(round(pct_repeats,0),"%")), size = 6, color = "black",
                position = position_stack(vjust = 0.5)) +
            theme_void() + scale_fill_hue(name="Repeats Location") +
            ggtitle("Deletions with homeology around repeats (both/one/no breakpoints)")
    ggsave(paste0(outDir,"/deletions_with_homeology.piechart.pdf"),width=5,height=5)

    dplot %>% fwrite(paste0(outDir, "/deletions_with_homeology.piechart.tsv"), sep="\t")

    # make barplot for deletions with homeology by repeat class
    dataDelHomeologyRepclass <- dataDelHomeologyRepclass %>% arrange(desc(prop))
    dataDelHomeologyRepclass$repClass <- factor(dataDelHomeologyRepclass$repClass,levels=as.character(unique(dataDelHomeologyRepclass$repClass)))

    p <- dataDelHomeologyRepclass %>% ggplot(aes(y=prop,x=repClass)) + 
        geom_bar(stat="identity",position=position_dodge(),fill='royalblue') + theme_prism(base_size=10) + 
        xlab('') + ylab('Proportion') + 
        ggtitle("Repeat class around deletions with homeology") + scale_fill_viridis_d(name="Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_discrete(drop=F)
    ggsave(paste0(outDir,"/deletions_with_homeology.repClass_barplot.pdf"),width=7,height=5)

    dataDelHomeologyRepclass %>% fwrite(paste0(outDir, "/deletions_with_homeology.repClass_barplot.tsv"), sep="\t")

}

make_piechart_for_proportion_of_repeat_regions_in_genome <- function(dataPath,outDir){

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    data <- fread(dataPath)

    dplot <- data.frame(
        group=c('Repeat region','Other region'),
        value=c(sum(data$size_prop)*100,100-sum(data$size_prop)*100)) %>%
        dplyr::mutate(group=factor(group,levels=c('Repeat region','Other region')))

    p <- dplot %>%
    ggplot(aes(x="",y=value,fill=group)) + geom_bar(stat="identity", width=1) + 
        coord_polar("y", start=0) + 
        geom_text(aes(label = paste0(round(value,0),"%")), size = 6, color = "black",
            position = position_stack(vjust = 0.5)) +
        theme_void() + scale_fill_hue(name="") +
        ggtitle("Repeat regions in the genome")

    ggsave(paste0(outDir,"/piechart_for_proportion_of_repeat_regions_in_genome.pdf"),p,width=5,height=5)

    dplot %>% fwrite(paste0(outDir, "/piechart_for_proportion_of_repeat_regions_in_genome.tsv"), sep="\t")

} 

make_boxplot_of_proportion_of_deletions_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2 <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena05_prop_of_dels_with_homeology_in_HRDetect_group_excluding_BRCA1_BRCA2")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # make boxplot
    dplot <- data %>%
        dplyr::filter(deletions != 0) %>%
        dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'), ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2'))

    p <- dplot %>%
        ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_boxplot() +
            theme_prism() +
            xlab("") +
            ylab("Proportion of deletions with homeology") +
            stat_compare_means(
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) + 
            coord_cartesian(ylim = c(0, 0.2))
        ggsave(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_HRDetect_high_low_excluding_BRCA1_BRCA2.ymax0.2.pdf"), width = 4, height = 5)

    p <- dplot %>%
        ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_boxplot() +
            theme_prism() +
            xlab("") +
            ylab("Proportion of deletions with homeology") +
            stat_compare_means(aes(label = ..p.format..),
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE),
                label.y=0.25
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5)) +
            coord_cartesian(ylim = c(0, 0.3))
        ggsave(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_HRDetect_high_low_excluding_BRCA1_BRCA2.pdf"), width = 4, height = 5)

    dplot %>% fwrite(paste0(outDir, "/boxplot_for_proportion_of_deletions_with_homeology_in_HRDetect_high_low_excluding_BRCA1_BRCA2.tsv"), sep="\t")

}


make_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2 <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena06_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2")
    ){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(dataPath) %>% 
        dplyr::filter(min_similarity==80, alignment_length_cutoff==30, ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2'))

    polq_exprs <- exprs %>% dplyr::filter(Name=='POLQ') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='POLQ_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(POLQ_nFPKM=as.numeric(as.character(POLQ_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(polq_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(POLQ_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/sample_data_with_POLQ_expression.tsv"),sep="\t")

    sample_info <- sample_info %>% dplyr::select(
        sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,POLQ_nFPKM,min_similarity,alignment_length_cutoff)

    ## make POLQ expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=POLQ_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("POLQ Gene Expression") + ggtitle("POLQ vs HRDetect") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/POLQ_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make POLQ boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_prism() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2.tsv"),sep="\t")


    ## make homeology rate vs POLQ plot  
    sample_info %>% fwrite(paste0(outDir,"/propDelHome_vs_POLQ_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Proportion of deletions with homeology") + ggtitle("propDelHome vs POLQ") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/propDelHome_vs_POLQ_expression.pdf"),p,width=5,height=5)

}


make_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2 <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena07_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2")
    ){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(dataPath) %>% 
        dplyr::filter(min_similarity==80, alignment_length_cutoff==30, ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2'))

    RAD52_exprs <- exprs %>% dplyr::filter(Name=='RAD52') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='RAD52_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(RAD52_nFPKM=as.numeric(as.character(RAD52_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(RAD52_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(RAD52_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/sample_data_with_RAD52_expression.tsv"),sep="\t")

    sample_info <- sample_info %>% dplyr::select(
        sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,RAD52_nFPKM,min_similarity,alignment_length_cutoff)

    ## make RAD52 expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=RAD52_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 vs HRDetect") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/RAD52_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make RAD52 boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_prism() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2.tsv"),sep="\t")


    ## make homeology rate vs RAD52 plot  
    sample_info %>% fwrite(paste0(outDir,"/propDelHome_vs_RAD52_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Proportion of deletions with homeology") + ggtitle("propDelHome vs RAD52") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/propDelHome_vs_RAD52_expression.pdf"),p,width=5,height=5)

}


make_POLQ_expression_in_BRCA2_and_control <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena08_POLQ_expression_in_BRCA2_and_control")
    ){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(dataPath) %>% 
        dplyr::filter(min_similarity==80, alignment_length_cutoff==30, allelic_status_brca1_brca2 %in% c('BRCA2','control'))

    polq_exprs <- exprs %>% dplyr::filter(Name=='POLQ') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='POLQ_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(POLQ_nFPKM=as.numeric(as.character(POLQ_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(polq_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(POLQ_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/sample_data_with_POLQ_expression.tsv"),sep="\t")

    sample_info <- sample_info %>% dplyr::select(
        sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,POLQ_nFPKM,min_similarity,alignment_length_cutoff)

    ## make POLQ expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=POLQ_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("POLQ Gene Expression") + ggtitle("POLQ vs HRDetect") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/POLQ_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make POLQ boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(deletions != 0)
    p <- dplot %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_prism() +
        stat_compare_means(comparisons=list(c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/POLQ_expression_in_BRCA2_and_control.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/POLQ_expression_in_BRCA2_and_control.tsv"),sep="\t")


    ## make homeology rate vs POLQ plot  
    sample_info %>% fwrite(paste0(outDir,"/propDelHome_vs_POLQ_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Proportion of deletions with homeology") + ggtitle("propDelHome vs POLQ") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/propDelHome_vs_POLQ_expression.pdf"),p,width=5,height=5)

}


make_RAD52_expression_in_BRCA2_and_control <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena09_RAD52_expression_in_BRCA2_and_control")
    ){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(dataPath) %>% 
        dplyr::filter(min_similarity==80, alignment_length_cutoff==30, allelic_status_brca1_brca2 %in% c('BRCA2','control'))

    RAD52_exprs <- exprs %>% dplyr::filter(Name=='RAD52') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='RAD52_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(RAD52_nFPKM=as.numeric(as.character(RAD52_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(RAD52_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(RAD52_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/sample_data_with_RAD52_expression.tsv"),sep="\t")

    sample_info <- sample_info %>% dplyr::select(
        sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,RAD52_nFPKM,min_similarity,alignment_length_cutoff)

    ## make RAD52 expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=RAD52_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression vs HRDetect") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/RAD52_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make RAD52 boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(deletions != 0)
    p <- dplot %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_prism() +
        stat_compare_means(comparisons=list(c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/RAD52_expression_in_BRCA2_and_control.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/RAD52_expression_in_BRCA2_and_control.tsv"),sep="\t")


    ## make homeology rate vs RAD52 plot  
    sample_info %>% fwrite(paste0(outDir,"/propDelHome_vs_RAD52_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Proportion of deletions with homeology") + ggtitle("propDelHome vs RAD52") + theme_prism() +
        stat_cor()
    ggsave(paste0(outDir,"/propDelHome_vs_RAD52_expression.pdf"),p,width=5,height=5)

}


make_boxplot_for_proportion_of_deletions_with_homeology_in_brca2_vs_control_for_different_deletion_lengths <- function(
    dataPath1KB_10KB=file.path(outdir,"/output_deletions_1KBto10KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath10KB_100KB=file.path(outdir,"/output_deletions_10KBto100KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath100KB_1MB=file.path(outdir,"/output_deletions_100KBto1MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath1MB_10MB=file.path(outdir,"/output_deletions_1MBto10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    dataPath10MB=file.path(outdir,"/output_deletions_above_10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena10_prop_of_dels_with_homeology_in_brca2_vs_control_for_different_deletion_lengths")){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data1KB_10KB <- fread(dataPath1KB_10KB) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30) %>% dplyr::mutate(deletion_size='1-10kb')
    data10KB_100KB <- fread(dataPath10KB_100KB) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30) %>% dplyr::mutate(deletion_size='10-100kb')
    data100KB_1MB <- fread(dataPath100KB_1MB) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30) %>% dplyr::mutate(deletion_size='100kb-1Mb')
    data1MB_10MB <- fread(dataPath1MB_10MB) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30) %>% dplyr::mutate(deletion_size='1Mb-10Mb')
    data10MB <- fread(dataPath10MB) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30) %>% dplyr::mutate(deletion_size='>10Mb')

    data <- rbind(data1KB_10KB,data10KB_100KB,data100KB_1MB,data1MB_10MB,data10MB) %>% 
        dplyr::mutate(deletion_size=factor(deletion_size,levels=c('1-10kb','10-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')))
    fwrite(data,paste0(outDir,"/sample_prop_of_dels_with_homeology_in_different_deletion_sizes.tsv"),sep="\t")

    # make boxplot
    dplot <- data %>% dplyr::filter(deletions != 0) %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))

    p <- dplot %>%
        ggplot(aes(x = allelic_status_brca1_brca2, y = homeology_rate)) +
            geom_boxplot() +
            theme_prism(base_size=12) +
            xlab("") +
            ylab("Proportion of deletions with homeology") +
            stat_compare_means(
                comparisons = list(c('BRCA2','control')),
                method = "wilcox.test", method.args = list("exact" = FALSE),
                label.y=0.32
            ) +
            coord_cartesian(ylim = c(0, 0.4)) +
            facet_wrap(~deletion_size,nrow=1)

    ggsave(paste0(outDir, "/boxplot_for_prop_of_dels_with_homeology_in_brca2_vs_control_for_different_deletion_lengths.ymax0.4.pdf"), width = 8, height = 3.5)

}


make_boxplots_for_BOPP_tumor_types <- function(
    dataPath=file.path(dataDir,"/BOPP_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    hrd_tbl_path="input/setton_hadi_choo_2023/hrd-supp-table_ZC.rds",
    outDir=file.path(figDir,"BOPP07_boxplots_for_BOPP_tumor_types")
){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    # Load data
    data <- fread(dataPath) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    # Load hrd table
    hrd_tbl <- readRDS(hrd_tbl_path)

    data <- data %>% dplyr::left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final),by=c('SAMPLE.TUMOR'='pair'))

    dplot <- data %>%
        dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'), tumor_type_final %in% c('BRCA','OV','PACA','PRAD'))

    p <- dplot %>% 
        ggplot(aes(x=tumor_type_final,y=homeology_rate)) + geom_boxplot() + theme_prism() + 
        coord_cartesian(ylim=c(0,0.5)) + xlab("") + ylab("Proportion of deletions with homeology")
    ggsave(paste0(outDir,"/boxplot_propDelHome_vs_tumor_type.pdf"),p,width=5,height=5)

    p <- dplot %>%
        ggplot(aes(x=tumor_type_final,y=homeology_rate)) + geom_boxplot() + theme_prism() + 
        coord_cartesian(ylim=c(0,0.5)) + xlab("") + ylab("Proportion of deletions with homeology") +
        facet_wrap(~allelic_status_brca1_brca2,nrow=1)
    ggsave(paste0(outDir,"/boxplot_propDelHome_vs_tumor_type.facet_BRCA2_control.pdf"),p,width=10,height=5)

    p <- dplot %>%
        ggplot(aes(x=allelic_status_brca1_brca2,y=homeology_rate)) + geom_boxplot() + theme_prism() + 
        coord_cartesian(ylim=c(0,0.5)) + 
        xlab("") + ylab("Proportion of deletions with homeology") +
        stat_compare_means(aes(label=..p.format..),method = "wilcox.test", method.args = list("exact" = FALSE), 
                label.y=0.45) +
        facet_wrap(~tumor_type_final,nrow=1)
    ggsave(paste0(outDir,"/boxplot_propDelHome_vs_tumor_type.facet_tumor_type.pdf"),p,width=8,height=5)

    dplot %>% fwrite(paste0(outDir,"/sample_propDelHome.tsv"),sep="\t")

    # calculate p-value
    dplot <- dplot %>% dplyr::mutate(group=paste0(allelic_status_brca1_brca2,' ',tumor_type_final))

    stat <- pairwise.wilcox.test(dplot$homeology_rate,dplot$group,p.adjust.method='none')$p.value %>% as.data.frame
    stat <- stat %>% rownames_to_column(var='group') %>% dplyr::select(group, everything())

    stat %>% fwrite(paste0(outDir,"/pairwise_wilcox_test.tsv"),sep="\t")



}



make_gene_expression_heatmap <- function(
    dataPath=file.path(dataDir,"/output_deletions_1KB-100KB_final/ssa_events_samples.tsv"),
    outDir=file.path(figDir,"Serena13_gene_expression_heatmap")
    ){

    # Make output directory
    if (!dir.exists(outDir)){
        dir.create(outDir, recursive = TRUE)
    }

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(dataPath) %>% 
        dplyr::filter(min_similarity==80, alignment_length_cutoff==30, allelic_status_brca1_brca2 %in% c('BRCA2','control'))

    
    polq_exprs <- exprs %>% dplyr::filter(Name=='POLQ') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='POLQ_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(POLQ_nFPKM=as.numeric(as.character(POLQ_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    RAD52_exprs <- exprs %>% dplyr::filter(Name=='RAD52') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='RAD52_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(RAD52_nFPKM=as.numeric(as.character(RAD52_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    brca2_exprs <- exprs %>% dplyr::filter(Name=='BRCA2') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='BRCA2_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(BRCA2_nFPKM=as.numeric(as.character(BRCA2_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% 
        left_join(polq_exprs %>% dplyr::rename('RNA_sampleid'='sampleid'),by=c('sampleid'='sample_name')) %>% # %>% dplyr::filter(!is.na(RAD52_nFPKM)) 
        left_join(RAD52_exprs %>% dplyr::select(-sampleid),by=c('sampleid'='sample_name')) %>% #%>% dplyr::filter(!is.na(RAD52_nFPKM))
        left_join(brca2_exprs %>% dplyr::select(-sampleid),by=c('sampleid'='sample_name')) %>% #%>% dplyr::filter(!is.na(BRCA2_nFPKM))
        dplyr::select(sampleid,RNA_sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,POLQ_nFPKM,RAD52_nFPKM,BRCA2_nFPKM,min_similarity,alignment_length_cutoff)

    sample_info %>% fwrite(paste0(outDir,"/sample_data_with_expression.tsv"),sep="\t")

    # make POLQ vs BRCA2 plot
    p <- sample_info %>% ggplot(aes(POLQ_nFPKM,BRCA2_nFPKM)) + geom_point() + stat_cor()
    ggsave(paste0(outDir,"/POLQ_vs_BRCA2.pdf"),p,width=5,height=5)
    p <- sample_info %>% ggplot(aes(POLQ_nFPKM,BRCA2_nFPKM)) + geom_point() + stat_cor() + facet_wrap(~allelic_status_brca1_brca2)
    ggsave(paste0(outDir,"/POLQ_vs_BRCA2.facet_genotype.pdf"),p,width=5,height=5)

    # prepare sample table
    sample_table <- sample_info %>% dplyr::filter(!is.na(POLQ_nFPKM),!is.na(RAD52_nFPKM),deletions>0)
    fwrite(sample_table,paste0(outDir,"/sample_table_for_heatmap.tsv"),sep="\t")

    # trim homeology rate
    print(summary(sample_table$homeology_rate))
    sample_table$homeology_rate <- ifelse(sample_table$homeology_rate>0.06,0.06,sample_table$homeology_rate)
    print(summary(sample_table$homeology_rate))

    # prepare gene experssion matrix
    exprs_mtx <- sample_table %>% dplyr::select(BRCA2_nFPKM,POLQ_nFPKM) %>% as.matrix()
    rownames(exprs_mtx) <- sample_table$sampleid
    exprs_mtx <- t(exprs_mtx)
    rownames(exprs_mtx) <- gsub('_nFPKM','',rownames(exprs_mtx))

    exprs_mtx <- t(scale(t(exprs_mtx)))

    exreme_value <- 2
    exprs_mtx[exprs_mtx>exreme_value] <- exreme_value
    exprs_mtx[exprs_mtx< -exreme_value] <- -exreme_value



    # make heatmap
    library(ComplexHeatmap)
    #set.seed(123)
    
    ## clustering by samples 
    library(circlize)
    col_fun = colorRamp2(c(0, 0.06), c("white", "purple"))

    column_ha = HeatmapAnnotation("Genotype"=sample_table$allelic_status_brca1_brca2,
                                  "HRDetect"=sample_table$HRDetect_group,
                                  "Proportion of deletions\nwith homeology"=sample_table$homeology_rate,
                                  col = list("Genotype" = c("BRCA2" = "#ff002b", "control" = "#0d0d0d"),
                                             "HRDetect" = c("HRDetect-high" = "#94e508", "HRDectect-intermediate"="#71a411", "HRDetect-low"="#498749"),
                                             "Proportion of deletions\nwith homeology" = col_fun
                                  ), annotation_name_gp = gpar(fontsize = 8)
                                )
    ht <- Heatmap(exprs_mtx,name="Expression level \n(Z-score)",show_column_names=FALSE,top_annotation = column_ha, cluster_rows = F)

    pdf(paste0(outDir,"/gene_expression_heatmap.clustered.pdf"),width=10,height=2)
    draw(ht)
    dev.off()

    ## sort by genotypes
    #set.seed(123)
    sample_table <- sample_table %>% dplyr::arrange(allelic_status_brca1_brca2,-HRDetect_score)
    exprs_mtx <- exprs_mtx[,sample_table$sampleid]

    library(circlize)
    col_fun = colorRamp2(c(0, 0.06), c("white", "purple"))

    column_ha = HeatmapAnnotation("Genotype"=sample_table$allelic_status_brca1_brca2,
                                  "HRDetect"=sample_table$HRDetect_group,
                                  "Proportion of deletions\nwith homeology"=sample_table$homeology_rate,
                                  col = list("Genotype" = c("BRCA2" = "#ff002b", "control" = "#0d0d0d"),
                                             "HRDetect" = c("HRDetect-high" = "#94e508", "HRDectect-intermediate"="#71a411", "HRDetect-low"="#498749"),
                                             "Proportion of deletions\nwith homeology" = col_fun
                                  ), annotation_name_gp = gpar(fontsize = 8)
                                )

    ht <- Heatmap(exprs_mtx,name="Expression level \n(Z-score)",show_column_names=FALSE,top_annotation = column_ha, cluster_rows = F, cluster_columns = F)

    pdf(paste0(outDir,"/gene_expression_heatmap.sort_by_genotypes.pdf"),width=10,height=2)
    draw(ht)
    dev.off()

    ## sort by HRDetect group
    #set.seed(123)
    print(table(sample_table$HRDetect_group))
    sample_table.X <- sample_table[sample(1:nrow(sample_table)),]
    sample_table.X <- sample_table.X %>% dplyr::mutate(HRDetect_group_order=ifelse(HRDetect_group=='HRDetect-high',1,ifelse(HRDetect_group=='HRDetect-low',3,2))) %>% 
        dplyr::arrange(HRDetect_group_order)
    exprs_mtx <- exprs_mtx[,sample_table.X$sampleid]

    sample_table.X %>% fwrite(paste0(outDir,"/sample_table_for_heatmap.sort_by_HRDetect_group.tsv"),sep="\t")

    library(circlize)
    col_fun = colorRamp2(c(0, 0.06), c("white", "purple"))

    column_ha = HeatmapAnnotation("Genotype"=sample_table.X$allelic_status_brca1_brca2,
                                  "HRDetect"=sample_table.X$HRDetect_group,
                                  "Proportion of deletions\nwith homeology"=sample_table.X$homeology_rate,
                                  col = list("Genotype" = c("BRCA2" = "#ff002b", "control" = "#0d0d0d"),
                                             "HRDetect" = c("HRDetect-high" = "#94e508", "HRDectect-intermediate"="#71a411", "HRDetect-low"="#498749"),
                                             "Proportion of deletions\nwith homeology" = col_fun
                                  ), annotation_name_gp = gpar(fontsize = 8)
                                )

    ht <- Heatmap(exprs_mtx,name="Expression level \n(Z-score)",show_column_names=FALSE,top_annotation = column_ha, cluster_rows = F, cluster_columns = F)

    pdf(paste0(outDir,"/gene_expression_heatmap.sort_by_HRDetect_group.pdf"),width=10,height=2)
    draw(ht)
    dev.off()


    ## sort by HRDetect group then by BRCA2 genotype
    #set.seed(123)
    sample_table.X <- sample_table %>% dplyr::mutate(HRDetect_group_order=ifelse(HRDetect_group=='HRDetect-high',1,ifelse(HRDetect_group=='HRDetect-low',3,2))) %>% 
        dplyr::arrange(HRDetect_group_order, allelic_status_brca1_brca2)
    exprs_mtx <- exprs_mtx[,sample_table.X$sampleid]

    sample_table.X %>% fwrite(paste0(outDir,"/sample_table_for_heatmap.sort_by_HRDetect_group_and_genotype.tsv"),sep="\t")

    library(circlize)
    col_fun = colorRamp2(c(0, 0.06), c("white", "purple"))

    column_ha = HeatmapAnnotation("Genotype"=sample_table.X$allelic_status_brca1_brca2,
                                  "HRDetect"=sample_table.X$HRDetect_group,
                                  "Proportion of deletions\nwith homeology"=sample_table.X$homeology_rate,
                                  col = list("Genotype" = c("BRCA2" = "#ff002b", "control" = "#0d0d0d"),
                                             "HRDetect" = c("HRDetect-high" = "#94e508", "HRDectect-intermediate"="#71a411", "HRDetect-low"="#498749"),
                                             "Proportion of deletions\nwith homeology" = col_fun
                                  ), annotation_name_gp = gpar(fontsize = 8)
                                )

    ht <- Heatmap(exprs_mtx,name="Expression level \n(Z-score)",show_column_names=FALSE,top_annotation = column_ha, cluster_rows = F, cluster_columns = F)

    pdf(paste0(outDir,"/gene_expression_heatmap.sort_by_HRDetect_group_and_genotype.pdf"),width=10,height=2)
    draw(ht)
    dev.off()


    ## sort by BRCA2
    #set.seed(123)
    sample_table <- sample_table %>% dplyr::arrange(BRCA2_nFPKM)
    exprs_mtx <- exprs_mtx[,sample_table$sampleid]

    library(circlize)
    col_fun = colorRamp2(c(0, 0.06), c("white", "purple"))

    column_ha = HeatmapAnnotation("Genotype"=sample_table$allelic_status_brca1_brca2,
                                  "HRDetect"=sample_table$HRDetect_group,
                                  "Proportion of deletions\nwith homeology"=sample_table$homeology_rate,
                                  col = list("Genotype" = c("BRCA2" = "#ff002b", "control" = "#0d0d0d"),
                                             "HRDetect" = c("HRDetect-high" = "#94e508", "HRDectect-intermediate"="#71a411", "HRDetect-low"="#498749"),
                                             "Proportion of deletions\nwith homeology" = col_fun
                                  ), annotation_name_gp = gpar(fontsize = 8)
                                )

    ht <- Heatmap(exprs_mtx,name="Expression level \n(Z-score)",show_column_names=FALSE,top_annotation = column_ha, cluster_rows = F, cluster_columns = F)

    pdf(paste0(outDir,"/gene_expression_heatmap.sort_by_BRCA2.pdf"),width=10,height=2)
    draw(ht)
    dev.off()


    

}
