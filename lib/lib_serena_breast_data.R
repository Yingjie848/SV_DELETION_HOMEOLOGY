examine_deletion_length.2 <- function(deletions, outDir,
        group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))){

        dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

        # histogram deletion length in BRCA1/BRCA2/control tumors
        p <- ggplot(deletions %>% dplyr::filter(get(group) %in% unlist(comparisons)), 
                aes(x=del_length)) + 
            geom_histogram(color="black", binwidth=0.2) + 
            scale_x_log10() +
            theme_bw() + xlab("Deletion length (bp)") + 
            facet_wrap( ~ get(group),scales="free_y")
        ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.pdf"),width=10.5,height=5)

        # boxplot deletion length in BRCA1/BRCA2/control tumors
        p <- ggplot(deletions %>% dplyr::filter(get(group) %in% unlist(comparisons)), 
                aes(x=get(group),y=del_length)) + 
            geom_violin() +
            geom_boxplot(outlier.shape = NA, width=0.1) + 
            geom_jitter(width=.25, size=0.3, alpha=0.1) + 
            scale_y_log10() +
            stat_compare_means(comparisons=comparisons) + 
            theme_bw() + xlab("") + ylab("Deletion length (bp)") 
        ggsave(paste0(outDir,"/deletion_len.boxplot_biallelic_vs_control.pdf"),width=3.5,height=5)


}


# make deletion-level plots for deletions with homeology
#analyze_ssa_homeology_events <- function(candidates,outDir){
make_deletion_level_plots <- function(candidates,outDir){

    dir.create(outDir,showWarnings = FALSE,recursive = TRUE)

    # homeology length histogram
    p <- candidates %>% ggplot(aes(x=alignment_length,y=after_stat(count)/(sum(after_stat(count))))) + 
            geom_histogram(color="black",binwidth=5) + ylab("Frequency") + xlab("Homeologous sequence length (bp)") + theme_bw()
    ggsave(paste0(outDir,"/homeology_len_hist.pdf"),width=7,height=5)

    # deletion length histogram
    p <- candidates %>% ggplot(aes(x=del_length,y=after_stat(count)/(sum(after_stat(count))))) + 
        geom_histogram(bins=50,color="black") + 
        scale_x_log10() + ylab("Frequency") + xlab("Deletion length (bp)") + theme_bw() 
        #facet_zoom(xlim = c(100, 50000), zoom.size = 1,split=TRUE)
    ggsave(paste0(outDir,"/deletion_len_hist.pdf"),width=7,height=5)


    # similarity histogram
    p <- candidates %>% ggplot(aes(x=identity,y=after_stat(count)/(sum(after_stat(count))))) + 
        geom_histogram(binwidth=1,color="black") + ylab("Frequency") + xlab("Similarity (%)") + theme_bw()
    ggsave(paste0(outDir,"/similarity_hist.pdf"),width=7,height=5)

    # scatterplot between deletion length and homeology length
    p <- candidates %>% ggplot(aes(x=alignment_length,y=del_length)) + geom_point() + stat_cor() + scale_y_log10() + xlab("Homeology length") + ylab("Deletion length") + theme_bw()
    ggsave(paste0(outDir,"/homeology_len_vs_deletion_len.pdf"),width=7,height=5)

    # scatterplot between deletion length and similarity
    p <- candidates %>% ggplot(aes(x=del_length,y=identity)) + geom_point() + scale_x_log10() + xlab("Deletion length") + ylab("Similarity") + theme_bw()
    ggsave(paste0(outDir,"/similarity_vs_deletion_len.pdf"),width=7,height=5)

    # scatterplot between homeology length and similarity
    p <- candidates %>% ggplot(aes(x=alignment_length,y=identity)) + geom_point() + stat_cor() + xlab("Homeology length") + ylab("Similarity") + theme_bw()
    ggsave(paste0(outDir,"/similarity_vs_homeology_len.pdf"),width=7,height=5)

}


# make deletion or sample-level plots for deletion length, alignment length, similarity for predicted deletions with homeology or avereaged by sample
# comparing between BRCA1/BRCA2/control and/or ER+/triple- tumors
# input: can be candidates or sample
make_deletion_or_sample_level_comparisons <- function(input,outDir, 
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))){

    dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

    input <- input %>% dplyr::filter(get(group) %in% unlist(comparisons))
    input %>% fwrite(paste0(outDir,"/candidates.tsv"),sep="\t")

    ## boxplot ---------------------------------------------------------------------------------------------------------------------------------

    # boxplot deletion length in BRCA1/BRCA2/control tumors
    p <- ggplot(input, aes(x=get(group),y=del_length)) + 
        geom_violin() +
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + scale_y_log10() + xlab("") + ylab("Deletion length (bp)") 
    ggsave(paste0(outDir,"/deletion_len.biallelic_vs_control.pdf"),width=3.5,height=5)

    # boxplot deletion length in BRCA1/BRCA2/control ER+/triple- tumors
    p <- ggplot(input, aes(x=get(group),y=del_length)) +
        geom_violin() + 
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) +  
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + scale_y_log10() + xlab("") + ylab("Deletion length (bp)") + 
        facet_wrap(~tumor_subtype)
    ggsave(paste0(outDir,"/deletion_len.biallelic_vs_control.ERpos_TripleNeg.pdf"),width=7,height=5)

    
    # boxplot homeology length in BRCA1/BRCA2/control tumors
    p <- ggplot(input, aes(x=get(group),y=alignment_length)) +
        geom_violin() + 
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + xlab("") + ylab("Homeologous sequence length (bp)") 
    ggsave(paste0(outDir,"/homeology_len.biallelic_vs_control.pdf"),width=3.5,height=5)

    # boxplot homeology length in BRCA1/BRCA2/control ER+/triple- tumors
    p <- ggplot(input, aes(x=get(group),y=alignment_length)) +
        geom_violin() + 
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + xlab("") + ylab("Homeologous sequence length (bp)") + 
        facet_wrap(~tumor_subtype)
    ggsave(paste0(outDir,"/homeology_len.biallelic_vs_control.ERpos_TripleNeg.pdf"),width=7,height=5)


    # boxplot similarity in BRCA1/BRCA2/control tumors
    p <- ggplot(input, aes(x=get(group),y=identity)) +
        geom_violin() + 
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + xlab("") + ylab("Similarity") 
    ggsave(paste0(outDir,"/similarity.biallelic_vs_control.pdf"),width=3.5,height=5)

    # boxplot similarity in BRCA1/BRCA2/control ER+/triple- tumors
    p <- ggplot(input, aes(x=get(group),y=identity)) +
        geom_violin() + 
        geom_boxplot(outlier.shape = NA, width=0.1) + geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        stat_compare_means(comparisons=comparisons) + 
        theme_bw() + xlab("") + ylab("Similarity") + 
        facet_wrap(~tumor_subtype)
    ggsave(paste0(outDir,"/similarity.biallelic_vs_control.ERpos_TripleNeg.pdf"),width=7,height=5)


    ## histograms ---------------------------------------------------------------------------------------------------------------------------------

    # histogram deletion length in BRCA1/BRCA2/control tumors
    p <- ggplot(input, 
            aes(x=del_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black", binwidth=0.2) + 
        theme_bw() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + ggtitle("Deletion length for deletions with homeology") + scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.pdf"),width=10.5,height=5)


    # normalized in each group
    p <- ggplot(input, 
            aes(x=del_length,y=stat(density*width))) +  # see https://stackoverflow.com/questions/22181132/normalizing-y-axis-in-histograms-in-r-ggplot-to-proportion-by-group
        geom_histogram(color="black", binwidth=0.2) + 
        theme_bw() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + ggtitle("Deletion length for deletions with homeology") + scale_x_log10() +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.normalizedInGroup.pdf"),width=10.5,height=5)


    # histogram deletion length in BRCA1/BRCA2/control tumors, group by ER+/triple-
    p <- ggplot(input %>% dplyr::mutate(group_name=paste0(tumor_subtype," / ",get(group))), 
            aes(x=del_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black", binwidth=0.2) + 
        theme_bw() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") + ggtitle("Deletion length for deletions with homeology") + scale_x_log10() +
        facet_wrap( ~ group_name,scales="fixed")
    ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.ERpos_TripleNeg.pdf"),width=10.5,height=5)

    # histogram homeology length in BRCA1/BRCA2/control tumors
    p <- ggplot(input, 
            aes(x=alignment_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black",binwidth=5) + 
        theme_bw() + xlab("Homeologous sequence length (bp)") + ylab("Proportion of deletions") + ggtitle("Homeologous sequence length for deletions with homeology") +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/homeology_len.histgram_biallelic_vs_control.pdf"),width=10.5,height=5)

    # normalized in each group
    p <- ggplot(input, 
            aes(x=alignment_length,y=after_stat(density*width))) + 
        geom_histogram(color="black",binwidth=5) + 
        theme_bw() + xlab("Homeologous sequence length (bp)") + ylab("Proportion of deletions") + ggtitle("Homeologous sequence length for deletions with homeology") +
        facet_wrap( ~ get(group),scales="fixed")
    ggsave(paste0(outDir,"/homeology_len.histgram_biallelic_vs_control.normalizedInGroup.pdf"),width=10.5,height=5)


    # histogram homeology length in BRCA1/BRCA2/control tumors, group by ER+/triple-
    p <- ggplot(input %>% dplyr::mutate(group_name=paste0(tumor_subtype," / ",get(group))), 
            aes(x=alignment_length,y=after_stat(count/sum(count)))) + 
        geom_histogram(color="black",binwidth=5) + 
        theme_bw() + xlab("Homeologous sequence length (bp)") + ylab("Proportion of deletions") + ggtitle("Homeologous sequence length for deletions with homeology") +
        facet_wrap( ~ group_name,scales="fixed")
    ggsave(paste0(outDir,"/homeology_len.histgram_biallelic_vs_control.ERpos_TripleNeg.pdf"),width=10.5,height=5)

    ## 

}


make_sample_level_comparisons_homeology_rate <- function(sample, outDir,
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))) {

        dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

        # make homeology_rate (n_events/deletions) boxplot 
        p <- sample %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_bialliec_BRCA_tumors.homeology_rate.pdf"), width = 3.5, height = 5)

        # make number of deletions with homeology boxplot 
        p <- sample %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = n_events)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Number of deletions with homeology") +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_bialliec_BRCA_tumors.deletions_with_homeology.pdf"), width = 3.5, height = 5)


        # make homeology_rate (n_events/deletions) boxplot by BRCA1/BRCA2/control
        # ymax = 0.5
        p <- sample %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE),
                label.y=0.4
            ) +
            coord_cartesian(ylim = c(0, 0.5)) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_bialliec_BRCA_tumors.homeology_rate.ymax0.5.pdf"), width = 3.5, height = 5)

        if(length(comparisons)==3){
            # make homeology_rate (n_events/deletions) boxplot by BRCA1/BRCA2/control
            # ymax = 0.5
            p <- sample %>%
                dplyr::filter(deletions != 0) %>%
                dplyr::filter(get(group) %in% unlist(comparisons)) %>%
                ggplot(aes(x = get(group), y = homeology_rate)) +
                geom_boxplot() +
                theme_bw() +
                xlab("") +
                ylab("Homeology rate") +
                stat_compare_means(
                    comparisons = comparisons,
                    method = "wilcox.test", method.args = list("exact" = FALSE),
                    label.y=c(0.3,0.36,0.42)
                ) +
                coord_cartesian(ylim = c(0, 0.5)) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
            ggsave(paste0(outDir, "/ssa_events_in_bialliec_BRCA_tumors.homeology_rate.ymax0.5.pdf"), width = 3.5, height = 5)
        }
        

        # make homeology rate violin
        p <- sample %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_bialliec_BRCA_tumors.homeology_rate.violin.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, hrd_status != "other") %>%
            ggplot(aes(x = hrd_status, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRD", "Non-HRD")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_hrd_tumors.homeology_rate.violin.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, hrd_status != "other", allelic_status_brca1_brca2 != 'BRCA1') %>%
            ggplot(aes(x = hrd_status, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRD", "Non-HRD")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_hrd_tumors_noBRCA1.homeology_rate.violin.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, HRDetect_group %in% c('HRDetect-high','HRDetect-low')) %>%
            ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
        ggsave(paste0(outDir, "/ssa_events_in_hrdetect_high_vs_low.homeology_rate.violin.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, HRDetect_group %in% c('HRDetect-high','HRDetect-low')) %>%
            ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
        ggsave(paste0(outDir, "/ssa_events_in_hrdetect_high_vs_low.homeology_rate.boxplot.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, HRDetect_group %in% c('HRDetect-high','HRDetect-low'), ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2')) %>%
            ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
        ggsave(paste0(outDir, "/ssa_events_in_hrdetect_high_vs_low_excluding_BRCA1_BRCA2.homeology_rate.violin.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, HRDetect_group %in% c('HRDetect-high','HRDetect-low'), ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2')) %>%
            ggplot(aes(x = HRDetect_group, y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRDetect-high", "HRDetect-low")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5))
        ggsave(paste0(outDir, "/ssa_events_in_hrdetect_high_vs_low_excluding_BRCA1_BRCA2.homeology_rate.boxplot.pdf"), width = 3.5, height = 5)


        p <- sample %>%
            dplyr::filter(deletions != 0, hrd_mut %in% c("Other HRD","Non-HRD")) %>%
            ggplot(aes(x = hrd_mut, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("Other HRD", "Non-HRD")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_other_hrd_tumors.homeology_rate.violin.pdf"), width = 3.5, height = 5)

        p <- sample %>%
            dplyr::filter(deletions != 0, hrd_status_hrdetect != "other") %>%
            ggplot(aes(x = hrd_status_hrdetect, y = homeology_rate)) +
            geom_violin(position = "identity") +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            stat_compare_means(
                comparisons = list(c("HRD", "Non-HRD")),
                method = "wilcox.test", method.args = list("exact" = FALSE)
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
        ggsave(paste0(outDir, "/ssa_events_in_hrd_tumors_hrdetect.homeology_rate.violin.pdf"), width = 3.5, height = 5)

}


# make plots for predicted SSA homeology candidates, compare between allelic status in different tumor subtype
make_sample_level_comparisons_other_details <- function(sample,outDir,
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))){

    dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

    # make plots
    make_deletion_or_sample_level_comparisons(input = sample,outDir=outDir,group=group,comparisons=comparisons)

}


# summarize candidates by sample
summarize_candidates <- function(candidates, full_deletions, outDir) {

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)


    # make sample level summary

    # get full tumor sample list, add no of events per sample, total number of deletions per sample, etc.
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
            tumor_subtype = tumor_subtype[1],
            allelic_status = allelic_status[1],
            allelic_status_brca1_brca2 = allelic_status_brca1_brca2[1],
            n_events = length(SAMPLE.TUMOR),
            del_length = mean(del_length),
            alignment_length = mean(alignment_length),
            identity = mean(identity)
        )
    sample_summary %>% fwrite(paste0(outDir, "/sample_summary.tsv"), sep = "\t")

    # get total number of deletions per sample
    deletions_per_sample <- full_deletions %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(deletions = length(SAMPLE.TUMOR))

    # merge sample_summary and deletions_per_sample for all samples
    sample <- get_full_tumor_list() %>%
        get_hrd() %>%
        merge(., sample_summary %>% dplyr::select(-tumor_subtype, -allelic_status, -allelic_status_brca1_brca2), by.x = "sampleid", by.y = "SAMPLE.TUMOR", all.x = TRUE) %>%
        merge(., deletions_per_sample, by.x = "sampleid", by.y = "SAMPLE.TUMOR", all.x = TRUE) %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events), deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        )
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

    # summarize samples by BRCA1/BRCA2/control, calculate proportion of homeology samples, proportion of homeology events
    sample %>%
        group_by(allelic_status_brca1_brca2) %>%
        summarise(total_samples = length(allelic_status_brca1_brca2), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors.tsv"), sep = "\t")

    # summarize samples by BRCA1/BRCA2/control, calculate proportion of homeology samples, proportion of homeology events
    sample_with_deletions <- sample %>% dplyr::filter(deletions > 0) 
    sample_with_deletions %>%
        group_by(allelic_status_brca1_brca2) %>%
        summarise(total_samples = length(allelic_status_brca1_brca2), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors_with_deletions.tsv"), sep = "\t")

    # add hrd_mut to candidates
    candidates <- candidates %>%
        merge(., sample %>% dplyr::select(sampleid, hrd_mut), by.x = "SAMPLE.TUMOR", by.y = "sampleid", all.x = TRUE)

    # make deletion level plots
    make_deletion_level_plots(candidates = candidates, outDir = paste0(outDir,"/deletion_level_plots"))

    # deletion-level comparisons for BRCA1/BRCA2/control and/or ER+/triple- tumors
    make_deletion_or_sample_level_comparisons(
        input=candidates %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')),
        outDir=paste0(outDir,"/deletion_level_comparisons"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

    # deletion-level comparisons for BRCA2/control and/or ER+/triple- tumors
    make_deletion_or_sample_level_comparisons(
        input=candidates %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control')),
        outDir=paste0(outDir,"/deletion_level_comparisons_BRCA2_control"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA2','control')))

    # deletion-level comparisons for BRCA2/Non-HRD and/or ER+/triple- tumors
    make_deletion_or_sample_level_comparisons(
        input=candidates %>% dplyr::filter(hrd_mut %in% c('BRCA2','Non-HRD')),
        outDir=paste0(outDir,"/deletion_level_comparisons_BRCA2_Non-HRD"),
        group="hrd_mut",
        comparisons=list(c('BRCA2','Non-HRD')))
    


    # make sample-level comparisons for homeology rate
    make_sample_level_comparisons_homeology_rate(
        sample = sample, 
        outDir = paste0(outDir, "/sample_level_comparisons_homeology_rate"),
        group = "allelic_status_brca1_brca2",
        comparisons = list(c("BRCA1", "BRCA2"), c("BRCA1", "control"), c("BRCA2", "control")))

    # make sample-level comparisons for homeology rate for BRCA2/control and ER+ tumors
    make_sample_level_comparisons_homeology_rate(
        sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2', 'control')), 
        outDir = paste0(outDir, "/sample_level_comparisons_BRCA2_control_homeology_rate"),
        group = "allelic_status_brca1_brca2",
        comparisons = list(c("BRCA2", "control")))

    # make sample-level comparisons for homeology rate for BRCA2/Non-HRD and ER+ tumors
    make_sample_level_comparisons_homeology_rate(
        sample = sample %>% dplyr::filter(hrd_mut %in% c('BRCA2', 'Non-HRD')), 
        outDir = paste0(outDir, "/sample_level_comparisons_BRCA2_Non-HRD_homeology_rate"),
        group = "hrd_mut",
        comparisons = list(c("BRCA2", "Non-HRD")))

    # make sample-level comparisons: averaged deletion length, homeology length, similarity for BRCA1/BRCA2/control and ER+/triple- tumors
    make_sample_level_comparisons_other_details(
        sample = sample,
        outDir = paste0(outDir,"/sample_level_comparisons_other_details"),
        group = "allelic_status_brca1_brca2",
        comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

    # make sample-level comparisons: averaged deletion length, homeology length, similarity for BRCA2/control and ER+/triple- tumors
    make_sample_level_comparisons_other_details(
        sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control')),
        outDir = paste0(outDir,"/sample_level_comparisons_BRCA2_control_other_details"),
        group = "allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control')))

    # make sample-level comparisons: averaged deletion length, homeology length, similarity for BRCA2/Non-HRD and ER+/triple- tumors
    make_sample_level_comparisons_other_details(
        sample = sample %>% dplyr::filter(hrd_mut %in% c('BRCA2','Non-HRD')),
        outDir = paste0(outDir,"/sample_level_comparisons_BRCA2_Non-HRD_other_details"),
        group = "hrd_mut",
        comparisons = list(c('BRCA2','Non-HRD')))

}


# get candidates for different similarities x alignment length, and then summarize candidates by samples
get_candidates <- function(full_deletions, blast_output_same_strand, outDir, 
    combineOnly = FALSE, 
    take_longest = TRUE, take_most_similar = FALSE, 
    similarities = c(80, 0, 90, 100), 
    alignment_lens = c(5, 10, 15, 20, 25, 30, 35, 40)) 
    {
    

    if (!combineOnly) {
        for (min_similarity in similarities) {
            cat("Similarity:", min_similarity, "\n")
            for (alignment_length_cutoff in alignment_lens) {
                cat("  Min alignment length:", alignment_length_cutoff, "\n")

                outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")
                dir.create(outdir_candidates, recursive = TRUE, showWarnings = FALSE)

                # for each event, get longest alignment with >=X% similarity
                candidates <- identify_ssa_candidates(blast_output_same_strand, min_similarity = min_similarity, alignment_length_cutoff = alignment_length_cutoff, take_longest = take_longest, take_most_similar = take_most_similar)
                # add deletion length
                candidates <- candidates %>% left_join(full_deletions %>% dplyr::select(SAMPLE.TUMOR, CHROM, start_position, end_position, del_length), by = c("SAMPLE.TUMOR", "CHROM", "start_position", "end_position"))
                # add tumor type and allelic status
                candidates <- add_tumor_type_and_allelic_status(candidates)
                candidates %>% fwrite(paste0(outdir_candidates, "/candidates.tsv"), sep = "\t")

                if (nrow(candidates) > 0) {
                    summarize_candidates(candidates = candidates, full_deletions = full_deletions, outDir = outdir_candidates)
                }
            }
        }
    }


    # combine candidates.tsv and ssa_events_for_each_sample.tsv files
    candidates <- data.frame()
    sample <- data.frame()
    for (min_similarity in similarities) {
        cat("Similarity:", min_similarity, "\n")
        for (alignment_length_cutoff in alignment_lens) {
            cat("  Min alignment length:", alignment_length_cutoff, "\n")

            outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")

            min_similarity_tmp <- min_similarity
            alignment_length_cutoff_tmp <- alignment_length_cutoff

            if (file.exists(paste0(outdir_candidates, "/ssa_events_for_each_sample.tsv"))) {
                candidates <- rbind(candidates, fread(paste0(outdir_candidates, "/candidates.tsv")) %>% 
                    dplyr::mutate(min_similarity = min_similarity_tmp, alignment_length_cutoff = alignment_length_cutoff_tmp))
                sample <- rbind(sample, fread(paste0(outdir_candidates, "/ssa_events_for_each_sample.tsv")) %>% 
                    dplyr::mutate(min_similarity = min_similarity_tmp, alignment_length_cutoff = alignment_length_cutoff_tmp))
            }
        }
    }
    list(candidates, sample)
}


# make sample-level plots by 2 factors: similarity and alignment length
make_sample_plots_similarity_vs_alignment_length <- function(
    sample,
    outDir="sample_plots_similarity_vs_alignment_length",
    similarities = c(0, 80, 90, 100), 
    alignment_lens = c(5, 10, 15, 20, 25, 30, 35, 40),
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))){

    dir.create(outDir)

    sample <- sample %>% dplyr::filter(min_similarity %in% similarities,alignment_length_cutoff %in% alignment_lens) %>%
        dplyr::mutate(min_similarity = factor(min_similarity),alignment_length_cutoff = factor(alignment_length_cutoff))

    # summarize by brca1/brca2/control groups
    summary <- sample %>% group_by(get(group),min_similarity,alignment_length_cutoff) %>% 
                                summarise(total_samples=length(get(group)),total_deletions=sum(deletions),homeology_samples=sum(n_events>=1,na.rm=TRUE),total_homeology_events=sum(n_events,na.rm=TRUE),mean_homeology_rate=mean(homeology_rate,na.rm=TRUE),median_homeology_rate=median(homeology_rate,na.rm=TRUE)) %>% 
                                dplyr::mutate(proportion_of_homeology_samples=homeology_samples/total_samples,
                                            mean_homeology_events_per_sample=total_homeology_events/total_samples,
                                            proportion_of_homeology_events=total_homeology_events/total_deletions
                                            )
    summary %>% fwrite(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.summary.tsv"),sep="\t")


    # homeology rate for min_similarity~alignment_length_cutoff in all samples
    p <- sample %>% dplyr::filter(deletions!=0,get(group) %in% unlist(comparisons)) %>%
                    ggplot(aes(x=get(group),y=homeology_rate)) + geom_boxplot() + theme_bw() + geom_jitter(size=0.3) + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        coord_cartesian(ylim=c(0,2)) +
                        facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.pdf"),width=12,height=10)


    p <- sample %>% dplyr::filter(deletions!=0,get(group) %in% unlist(comparisons)) %>%
                    ggplot(aes(x=get(group),y=homeology_rate)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        coord_cartesian(ylim=c(0,1.5)) +
                        facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.violin.pdf"),width=12,height=10)


    p <- sample %>% dplyr::filter(deletions!=0,hrd_status != "other") %>%
                    ggplot(aes(x=hrd_status,y=homeology_rate)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=list(c('HRD','Non-HRD')), 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        coord_cartesian(ylim=c(0,1.5)) +
                        facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_hrd_tumors.homeology_rate.violin.pdf"),width=12,height=10)


    p <- sample %>% dplyr::filter(deletions!=0,hrd_mut %in% c("Other HRD","Non-HRD")) %>%
                    dplyr::mutate(hrd_mut = factor(hrd_mut,levels=c("Other HRD","Non-HRD"))) %>%
                    ggplot(aes(x=hrd_mut,y=homeology_rate)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=list(c('Other HRD','Non-HRD')), 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        coord_cartesian(ylim=c(0,1.5)) +
                        facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_other_hrd_tumors.homeology_rate.violin.pdf"),width=12,height=10)


    if('HRDetect_score' %in% colnames(sample)){
        # HRDetect score
        p <- sample %>% dplyr::filter(deletions!=0,get(group) %in% unlist(comparisons)) %>%
                        ggplot(aes(x=get(group),y=HRDetect_score)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("HRDetect score") + ggtitle("HRDetect score") + 
                            stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                                method="wilcox.test",method.args=list('exact'=FALSE)) +
                            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                            coord_cartesian(ylim=c(0,1.5)) +
                            facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score.violin.pdf"),width=12,height=10)


        # homeology_rate vs HRDetect score in brca1 samples
        p <- sample %>% dplyr::filter(deletions!=0,get(group)=="BRCA1") %>%
                        ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (BRCA1)") +
                            facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_brca1.pdf"),width=12,height=10)

        # homeology_rate vs HRDetect score in brca2 samples
        p <- sample %>% dplyr::filter(deletions!=0,get(group)=="BRCA2") %>%
                        ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (BRCA2)") +
                            facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_brca2.pdf"),width=12,height=10)

        # homeology_rate vs HRDetect score in control samples
        p <- sample %>% dplyr::filter(deletions!=0,get(group)=="control") %>%
                        ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (control)") +
                            facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_control.pdf"),width=12,height=10)

    }
}


# 
make_sample_plots_alignment_length <- function(sample,outDir,
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))
    {

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    # summarize by brca1/brca2/control groups
    summary <- sample %>% group_by(get(group),min_similarity,alignment_length_cutoff) %>%
                                summarise(total_samples=length(get(group)),total_deletions=sum(deletions),homeology_samples=sum(n_events>=1,na.rm=TRUE),total_homeology_events=sum(n_events,na.rm=TRUE),mean_homeology_rate=mean(homeology_rate,na.rm=TRUE),median_homeology_rate=median(homeology_rate,na.rm=TRUE)) %>% 
                                dplyr::mutate(proportion_of_homeology_samples=homeology_samples/total_samples,
                                            mean_homeology_events_per_sample=total_homeology_events/total_samples,
                                            proportion_of_homeology_events=total_homeology_events/total_deletions
                                            )
    summary %>% fwrite(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.summary.tsv"),sep="\t")


    p <- sample %>% dplyr::filter(deletions!=0, get(group) %in% unlist(comparisons)) %>%
                    ggplot(aes(x=get(group),y=homeology_rate)) + geom_boxplot(outlier.shape=NA) + 
                        theme_bw() + #geom_jitter(size=0.3, width=0.5) + 
                        xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        facet_grid(~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.pdf"),width=7,height=5)


    p <- sample %>% dplyr::filter(deletions!=0, get(group) %in% unlist(comparisons)) %>%
                    ggplot(aes(x=get(group),y=homeology_rate)) + 
                        geom_violin(position = "identity") + 
                        geom_sina(size=0.3) + 
                        theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        #coord_cartesian(ylim=c(0,1.5)) +
                        facet_wrap(~alignment_length_cutoff,scales="free_y",nrow=1) 
    ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.violin.pdf"),width=10,height=5)


    p <- sample %>% dplyr::filter(deletions!=0, get(group) %in% unlist(comparisons)) %>%
                    ggplot(aes(x=get(group),y=deletions)) + 
                        geom_violin(position = "identity") + 
                        geom_sina(size=0.3) + 
                        theme_bw() + xlab("") + ylab("Deletions") + ggtitle("Number of deletions") + 
                        stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        #coord_cartesian(ylim=c(0,1.5)) +
                        facet_wrap(~alignment_length_cutoff,scales="free_y",nrow=1) 
    ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.deletions.violin.pdf"),width=10,height=5)


    p <- sample %>% dplyr::filter(deletions!=0,hrd_status != "other") %>%
                    ggplot(aes(x=hrd_status,y=homeology_rate)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=list(c('HRD','Non-HRD')), 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        #coord_cartesian(ylim=c(0,1.5)) +
                        facet_grid(~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_hrd_tumors.homeology_rate.violin.pdf"),width=7,height=5)


    p <- sample %>% dplyr::filter(deletions!=0,hrd_mut %in% c("Other HRD","Non-HRD")) %>%
                    dplyr::mutate(hrd_mut = factor(hrd_mut,levels=c("Other HRD","Non-HRD"))) %>%
                    ggplot(aes(x=hrd_mut,y=homeology_rate)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                        stat_compare_means(aes(label=p.signif),comparisons=list(c('Other HRD','Non-HRD')), 
                                            method="wilcox.test",method.args=list('exact'=FALSE)) +
                        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                        #coord_cartesian(ylim=c(0,1.5)) +
                        facet_grid(~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/ssa_events_in_other_hrd_tumors.homeology_rate.violin.pdf"),width=7,height=5)


    if('HRDetect_score' %in% colnames(sample)){
        # HRDetect score
        p <- sample %>% dplyr::filter(deletions!=0, get(group) %in% unlist(comparisons)) %>%
                        ggplot(aes(x=get(group),y=HRDetect_score)) + geom_violin(position = "identity") + geom_sina(size=0.3) + theme_bw() + xlab("") + ylab("HRDetect score") + ggtitle("HRDetect score") + 
                            stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                                method="wilcox.test",method.args=list('exact'=FALSE)) +
                            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                            #coord_cartesian(ylim=c(0,1.5)) +
                            facet_grid(~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score.violin.pdf"),width=7,height=5)


        p <- sample %>% dplyr::filter(deletions!=0, HRDetect_group %in% c('HRDetect-high','HRDetect-low')) %>%
                        ggplot(aes(x=HRDetect_group,y=homeology_rate)) + geom_violin(position = "identity") + 
                            geom_sina(size=0.3) + theme_bw() + 
                            xlab("") + ylab("Homeology rate") + ggtitle("Homeology rate") + 
                            stat_compare_means(aes(label=p.signif),comparisons=comparisons, 
                                                method="wilcox.test",method.args=list('exact'=FALSE)) +
                            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1)) +
                            facet_grid(~alignment_length_cutoff,scales="free") 
        ggsave(paste0(outDir,"/ssa_events_in_hrdetect_high_vs_low_tumors.homeology_rate.violin.pdf"),width=7,height=5)


        # homeology_rate vs HRDetect score in brca1 samples
        if('BRCA1' %in% sample[,get(group)]){
            p <- sample %>% dplyr::filter(deletions!=0,get(group)=="BRCA1") %>%
                            ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (BRCA1)") +
                                facet_grid(~alignment_length_cutoff,scales="free") 
            ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_brca1.pdf"),width=7,height=5)
        }

        # homeology_rate vs HRDetect score in brca2 samples
        if('BRCA2' %in% sample[,get(group)]){
            p <- sample %>% dplyr::filter(deletions!=0,get(group)=="BRCA2") %>%
                            ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (BRCA2)") +
                                facet_grid(~alignment_length_cutoff,scales="free") 
            ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_brca2.pdf"),width=7,height=5)
        }

        # homeology_rate vs HRDetect score in control samples
        if('control' %in% sample[,get(group)]){
            p <- sample %>% dplyr::filter(deletions!=0,get(group)=="control") %>%
                            ggplot(aes(x=HRDetect_score,y=homeology_rate)) + geom_point() + theme_bw() + xlab("HRDetect score") + ylab("Homeology rate") + ggtitle("Homeology rate vs HRDetect score (control)") +
                                facet_grid(~alignment_length_cutoff,scales="free") 
            ggsave(paste0(outDir,"/ssa_events_in_bialliec_brca1_brca2_tumors.hrdetect_score_vs_homeology_rate_in_control.pdf"),width=7,height=5)
        }

    }
}


count_number_of_dels_homeology <- function(outDir) {
    blast_summary <- fread(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors.summary.tsv")) %>%
        dplyr::filter(allelic_status_brca1_brca2 != "other")

    summary <- blast_summary %>%
        dplyr::filter(alignment_length_cutoff != 1) %>%
        group_by(min_similarity, alignment_length_cutoff) %>%
        summarise(total_samples = sum(total_samples), total_deletions = sum(total_deletions), total_homeology_events = sum(total_homeology_events), homeology_rate = total_homeology_events / total_deletions * 100) %>%
        dplyr::arrange(min_similarity, alignment_length_cutoff) %>%
        dplyr::rename(minimal_homeology_length = alignment_length_cutoff)

    summary %>% fwrite(paste0(outDir, "/count_number_of_dels_homeology.tsv"), sep = "\t")
}


# get candidates for short deleitons in similarities x alignment length, and then summarize candidates by samples
get_candidates_for_short_deletions <- function(full_deletions, blast_output_same_strand, similarities = c(0, 80, 90, 100), alignment_lens = 4:9, outDir, combineOnly = FALSE) {
    if (!combineOnly) {
        for (min_similarity in similarities) {
            cat("Similarity:", min_similarity, "\n")
            for (alignment_length_cutoff in alignment_lens) {
                cat("  Min alignment length:", alignment_length_cutoff, "\n")

                outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")
                dir.create(outdir_candidates)

                # for each event, get longest alignment with >=X% similarity
                candidates <- identify_ssa_candidates(blast_output_same_strand, min_similarity = min_similarity, alignment_length_cutoff = alignment_length_cutoff)
                # add deletion length
                candidates <- candidates %>% left_join(full_deletions %>% dplyr::select(SAMPLE.TUMOR, CHROM, start_position, end_position, del_length), by = c("SAMPLE.TUMOR", "CHROM", "start_position", "end_position"))
                # add tumor type and allelic status
                candidates <- add_tumor_type_and_allelic_status(candidates)
                candidates %>% fwrite(paste0(outdir_candidates, "/candidates.tsv"), sep = "\t")

                if (nrow(candidates) > 0) {
                    summarize_candidates(candidates = candidates, full_deletions = full_deletions, outDir = outdir_candidates)
                }
            }
        }
    }


    # combine candidates.tsv files and compare BRCA1/BRCA2/control
    candidates <- data.frame()
    sample <- data.frame()
    for (min_similarity in similarities) {
        cat("Similarity:", min_similarity, "\n")
        for (alignment_length_cutoff in alignment_lens) {
            cat("  Min alignment length:", alignment_length_cutoff, "\n")

            outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")

            min_similarity_tmp <- min_similarity
            alignment_length_cutoff_tmp <- alignment_length_cutoff

            if (file.exists(paste0(outdir_candidates, "/ssa_events_for_each_sample.tsv"))) {
                candidates <- rbind(candidates, fread(paste0(outdir_candidates, "/candidates.tsv")) %>% dplyr::mutate(min_similarity = min_similarity_tmp, alignment_length_cutoff = alignment_length_cutoff_tmp))
                sample <- rbind(sample, fread(paste0(outdir_candidates, "/ssa_events_for_each_sample.tsv")) %>% dplyr::mutate(min_similarity = min_similarity_tmp, alignment_length_cutoff = alignment_length_cutoff_tmp))
            }
        }
    }
    list(candidates, sample)
}


# compare alignment length for most similar and longest alignment
compare_alignment_length_most_similar_and_longest <- function(candidates,outDir){
    candidates_most_similar <- candidates %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length,alignment_length,identity,min_similarity,alignment_length_cutoff,allelic_status_brca1_brca2) %>%
                                        dplyr::rename(alignment_length_most_similar=alignment_length,similarity_most_similar=identity)
    candidates_longest <- fread(paste0(outdir,"/output_deletions_above_1KB/ssa_events_candidates.tsv")) %>% 
                        dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length,alignment_length,identity,min_similarity,alignment_length_cutoff,allelic_status_brca1_brca2) %>% 
                        dplyr::rename(alignment_length_longest=alignment_length,similarity_longest=identity)
    candidates_comb <- candidates_most_similar %>% left_join(candidates_longest) %>% dplyr::filter(alignment_length_cutoff!=1)

    p <- candidates_comb %>% ggplot(aes(x=alignment_length_longest,y=alignment_length_most_similar)) + geom_point() + theme_bw() + 
                                ggtitle("Homeology sequence length (most similar vs longest)") + stat_cor(aes(label = ..r.label..)) +
                                facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/most_similar_vs_longest.alignment_length.pdf"),width=12,height=10)

    p <- candidates_comb %>% ggplot(aes(x=similarity_longest,y=similarity_most_similar)) + geom_point() + theme_bw() + 
                                ggtitle("Similarity (most similar vs longest)") + stat_cor(aes(label = ..r.label..)) +
                                facet_grid(min_similarity~alignment_length_cutoff,scales="free") 
    ggsave(paste0(outDir,"/most_similar_vs_longest.similarity.pdf"),width=12,height=10)
}

# check POLQ expression in HRDetect high and low groups
check_POLQ_expression_in_HRDetect_group <- function(outDir=outDir){

    dir.create(paste0(outDir,"/plots_POLQ_expression"))

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(paste0(outDir,"/ssa_events_samples.tsv")) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    polq_exprs <- exprs %>% dplyr::filter(Name=='POLQ') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='POLQ_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(POLQ_nFPKM=as.numeric(as.character(POLQ_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(polq_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(POLQ_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/plots_POLQ_expression/sample_data_with_POLQ_expression.tsv"),sep="\t")
    sample_info <- sample_info %>% dplyr::select(sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,POLQ_nFPKM,min_similarity,alignment_length_cutoff)

    ## make POLQ expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=POLQ_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression vs Homeology rate") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make POLQ boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_HRDetect_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_HRDetect_group.tsv"),sep="\t")

    # BRCA1/BRCA2/control group
    dplot <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control'))
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_BRCA1_BRCA2_genotype_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_BRCA1_BRCA2_genotype_group.tsv"),sep="\t")

    # BRCA2/control group
    dplot <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control')) %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_BRCA2_genotype_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_POLQ_expression/POLQ_expression_in_BRCA2_genotype_group.tsv"),sep="\t")


    ## make homeology rate vs POLQ plot 
    sample_info %>% fwrite(paste0(outDir,"/plots_POLQ_expression/homeology_rate_vs_POLQ_expression.tsv"),sep="\t")

    # 
    p <- sample_info %>%
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs POLQ Gene Expression") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_POLQ_expression/homeology_rate_vs_POLQ_expression.pdf"),p,width=5,height=5)
    
    # HRDetect group
    p <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low')) %>% 
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs POLQ Gene Expression") + theme_classic2() +
        stat_cor() +
        facet_wrap(~HRDetect_group,scales="free")
    ggsave(paste0(outDir,"/plots_POLQ_expression/homeology_rate_vs_POLQ_expression_in_HRDetect_group.pdf"),p,width=10,height=5)

    # BRCA1/BRCA2/control group
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>% 
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs POLQ Gene Expression") + theme_classic2() +
        stat_cor() +
        facet_wrap(~allelic_status_brca1_brca2,scales="free")
    ggsave(paste0(outDir,"/plots_POLQ_expression/homeology_rate_vs_POLQ_expression_in_BRCA1_BRCA2_genotype_group.pdf"),p,width=10,height=5)


}


# check POLQ expression in HRDetect high and low groups excluding BRCA1 and BRCA2
check_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2 <- function(outDir=outDir){

    dir.create(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2"))

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(paste0(outDir,"/ssa_events_samples.tsv")) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30, ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2'))

    polq_exprs <- exprs %>% dplyr::filter(Name=='POLQ') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='POLQ_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(POLQ_nFPKM=as.numeric(as.character(POLQ_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(polq_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(POLQ_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/sample_data_with_POLQ_expression.tsv"),sep="\t")
    sample_info <- sample_info %>% dplyr::select(sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,POLQ_nFPKM,min_similarity,alignment_length_cutoff)

    ## make POLQ expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=POLQ_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression vs Homeology rate") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/POLQ_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make POLQ boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=POLQ_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("POLQ Gene Expression") + ggtitle("POLQ Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/POLQ_expression_in_HRDetect_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/POLQ_expression_in_HRDetect_group.tsv"),sep="\t")


    ## make homeology rate vs POLQ plot  
    sample_info %>% fwrite(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/homeology_rate_vs_POLQ_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=POLQ_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("POLQ Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs POLQ Gene Expression") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_POLQ_expression_excluding_BRCA1_BRCA2/homeology_rate_vs_POLQ_expression.pdf"),p,width=5,height=5)

}


# check RAD52 expression in HRDetect high and low groups
check_RAD52_expression_in_HRDetect_group <- function(outDir=outDir){

    dir.create(paste0(outDir,"/plots_RAD52_expression"))

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(paste0(outDir,"/ssa_events_samples.tsv")) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30)

    RAD52_exprs <- exprs %>% dplyr::filter(Name=='RAD52') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='RAD52_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(RAD52_nFPKM=as.numeric(as.character(RAD52_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(RAD52_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(RAD52_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/plots_RAD52_expression/sample_data_with_RAD52_expression.tsv"),sep="\t")
    sample_info <- sample_info %>% dplyr::select(sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,RAD52_nFPKM,min_similarity,alignment_length_cutoff)

    ## make RAD52 expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=RAD52_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression vs Homeology rate") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make RAD52 boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_HRDetect_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_HRDetect_group.tsv"),sep="\t")

    # BRCA1/BRCA2/control group
    dplot <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control'))
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_BRCA1_BRCA2_genotype_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_BRCA1_BRCA2_genotype_group.tsv"),sep="\t")

    # BRCA2/control group
    dplot <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control')) %>% 
        ggplot(aes(x=allelic_status_brca1_brca2,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('BRCA2','control')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_BRCA2_genotype_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_RAD52_expression/RAD52_expression_in_BRCA2_genotype_group.tsv"),sep="\t")


    ## make homeology rate vs RAD52 plot 
    sample_info %>% fwrite(paste0(outDir,"/plots_RAD52_expression/homeology_rate_vs_RAD52_expression.tsv"),sep="\t")

    # 
    p <- sample_info %>%
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs RAD52 Gene Expression") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_RAD52_expression/homeology_rate_vs_RAD52_expression.pdf"),p,width=5,height=5)
    
    # HRDetect group
    p <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low')) %>% 
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs RAD52 Gene Expression") + theme_classic2() +
        stat_cor() +
        facet_wrap(~HRDetect_group,scales="free")
    ggsave(paste0(outDir,"/plots_RAD52_expression/homeology_rate_vs_RAD52_expression_in_HRDetect_group.pdf"),p,width=10,height=5)

    # BRCA1/BRCA2/control group
    p <- sample_info %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>% 
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs RAD52 Gene Expression") + theme_classic2() +
        stat_cor() +
        facet_wrap(~allelic_status_brca1_brca2,scales="free")
    ggsave(paste0(outDir,"/plots_RAD52_expression/homeology_rate_vs_RAD52_expression_in_BRCA1_BRCA2_genotype_group.pdf"),p,width=10,height=5)


}


# check RAD52 expression in HRDetect high and low groups excluding BRCA1 and BRCA2
check_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2 <- function(outDir=outDir){

    dir.create(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2"))

    exprs <- readRDS("input/Serena_RNAseq_fpkm_table.rds")
    exprs_sample_info <- fread("input/Serena_data/sample_info.csv")
    sample_info <- fread(paste0(outDir,"/ssa_events_samples.tsv")) %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30, ! allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2'))

    RAD52_exprs <- exprs %>% dplyr::filter(Name=='RAD52') %>% dplyr::select(-UNIQID,-Ensembl,-Name,-Source,-loc) %>%
        pivot_longer(cols=colnames(exprs)[6:ncol(exprs)],names_to='sampleid',values_to='RAD52_nFPKM') %>%
        left_join(exprs_sample_info %>% dplyr::select(RNA.name,sample_name),by=c('sampleid'='RNA.name')) %>%
        dplyr::mutate(RAD52_nFPKM=as.numeric(as.character(RAD52_nFPKM))) %>%
        dplyr::filter(!is.na(sample_name)) %>% droplevels()

    sample_info <- sample_info %>% left_join(RAD52_exprs,by=c('sampleid'='sample_name')) %>% dplyr::filter(!is.na(RAD52_nFPKM))

    sample_info %>% fwrite(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/sample_data_with_RAD52_expression.tsv"),sep="\t")
    sample_info <- sample_info %>% dplyr::select(sampleid,allelic_status_brca1_brca2,HRDetect_score,HRDetect_group,deletions,n_events,homeology_rate,RAD52_nFPKM,min_similarity,alignment_length_cutoff)

    ## make RAD52 expression vs HRDetect score scatterplot
    p <- sample_info %>%
        ggplot(aes(x=HRDetect_score,y=RAD52_nFPKM)) + 
        geom_point() + xlab("HRDetect score") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression vs Homeology rate") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/RAD52_expression_vs_HRDetectScore.pdf"),p,width=5,height=5)
    

    ## make RAD52 boxplot

    # HRDetect group
    dplot <- sample_info %>% dplyr::filter(HRDetect_group %in% c('HRDetect-high','HRDetect-low'))
    p <- dplot %>% 
        ggplot(aes(x=HRDetect_group,y=RAD52_nFPKM)) + 
        geom_boxplot() + xlab("") + ylab("RAD52 Gene Expression") + ggtitle("RAD52 Gene Expression") + theme_classic2() +
        stat_compare_means(comparisons=list(c('HRDetect-high','HRDetect-low')),method='wilcox.test') 

    ggsave(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/RAD52_expression_in_HRDetect_group.pdf"),p,width=4,height=5)
    dplot %>% fwrite(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/RAD52_expression_in_HRDetect_group.tsv"),sep="\t")


    ## make homeology rate vs RAD52 plot  
    sample_info %>% fwrite(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/homeology_rate_vs_RAD52_expression.tsv"),sep="\t")
    p <- sample_info %>%
        ggplot(aes(x=RAD52_nFPKM,y=homeology_rate)) + 
        geom_point() + xlab("RAD52 Gene Expression") + ylab("Homeology rate") + ggtitle("Homeology rate vs RAD52 Gene Expression") + theme_classic2() +
        stat_cor()
    ggsave(paste0(outDir,"/plots_RAD52_expression_excluding_BRCA1_BRCA2/homeology_rate_vs_RAD52_expression.pdf"),p,width=5,height=5)

}
