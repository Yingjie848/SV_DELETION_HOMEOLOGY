# make deletion-level plots for deletions with homeology, including homeology length histogram, deletion length histogram, similarity histogram, scatterplot between deletion length and homeology length, scatterplot between deletion length and similarity, scatterplot between homeology length and similarity

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



count_number_of_dels_homeology <- function(){

    blast_summary <- fread("output/serena_SV_deletions_marcin_code/ssa_events_in_bialliec_brca1_brca2_tumors.summary.tsv") %>% 
                     dplyr::filter(allelic_status_brca1_brca2!='other')
    
    summary <- blast_summary %>% dplyr::filter(alignment_length_cutoff!=1) %>%
                group_by(min_similarity,alignment_length_cutoff) %>% 
                summarise(total_samples=sum(total_samples),total_deletions=sum(total_deletions),total_homeology_events=sum(total_homeology_events),homeology_rate=total_homeology_events/total_deletions*100) %>%
                dplyr::arrange(min_similarity,alignment_length_cutoff)

    summary %>% fwrite("output/serena_SV_deletions_marcin_code/count_number_of_dels_homeology.tsv",sep="\t")

} 


# check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
make_distance_to_breakpoint_plot <- function(candidates,outDir,
    similarities = c(0,80,90,100),
    alignment_lens = c(5,10,15,20,25,30,35,40)){


    for(thres_min_similarity in similarities){
        cat("Similarity:",thres_min_similarity,"\n")
        for(thres_alignment_length_cutoff in alignment_lens){
            cat("  Min alignment length:",thres_alignment_length_cutoff,"\n")

            outdir_candidates <- paste0(outDir, "/results_similarity_", thres_min_similarity, "pct_", thres_alignment_length_cutoff, "bp")
            dir.create(outdir_candidates,showWarnings = FALSE,recursive = TRUE)

            candidates_tmp <- candidates %>% 
                dplyr::mutate(q_start_dbrk=q_start-100,q_end_dbrk=q_end-100,s_start_dbrk=q_end-100,s_end_dbrk=s_end-100) %>%
                dplyr::mutate(
                    q_dbrk=ifelse(abs(q_start_dbrk)<abs(q_end_dbrk),q_start_dbrk,q_end_dbrk),
                    s_dbrk=ifelse(abs(s_start_dbrk)<abs(s_end_dbrk),s_start_dbrk,s_end_dbrk),
                    dbrk=ifelse(q_dbrk<s_dbrk, q_dbrk, s_dbrk)
                ) %>%
                dplyr::filter(min_similarity==thres_min_similarity, alignment_length_cutoff==thres_alignment_length_cutoff)
            print(nrow(candidates_tmp))
            candidates_tmp %>% fwrite(paste0(outdir_candidates,"/candidates.add_dbrk.tsv"),sep="\t")

            p1 <- candidates_tmp %>%
                ggplot(aes(x=q_dbrk,y=after_stat(count/sum(count)))) + geom_histogram(color="black",binwidth=25) + xlab("Distance to breakpoint") + ylab("Frequency") + ggtitle("Distance to 1st breakpoint") + theme_bw()
            p2 <- candidates_tmp %>% 
                ggplot(aes(x=s_dbrk,y=after_stat(count/sum(count)))) + geom_histogram(color="black",binwidth=25) + xlab("Distance to breakpoint") + ylab("Frequency") + ggtitle("Distance to 2nd breakpoint") + theme_bw()
            p3 <- candidates_tmp %>% 
                ggplot(aes(x=dbrk,y=after_stat(count/sum(count)))) + geom_histogram(color="black",binwidth=25) + xlab("Distance to breakpoint") + ylab("Frequency") + ggtitle("Distance to breakpoints") + theme_bw()

            p <- ggarrange(p1,p2,p3,nrow=1)
            ggsave(paste0(outdir_candidates,"/distance_to_breakpoints.pdf"),width=10.5,height=3.5)

        }
    }
}


examine_deletion_length <- function(deletions, outDir,
        group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control'))){

        dir.create(outDir,recursive=TRUE,showWarnings=FALSE)

        # histogram deletion length in BRCA1/BRCA2/control tumors
        p <- ggplot(deletions %>% dplyr::filter(get(group) %in% unlist(comparisons)), 
                aes(x=del_length,y=after_stat(count/sum(count)))) + 
            geom_histogram(color="black", binwidth=0.2) + 
            scale_x_log10() +
            theme_bw() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") +
            facet_wrap( ~ get(group),scales="fixed")
        ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.pdf"),width=10.5,height=5)

        # normalize in each group
        p <- ggplot(deletions %>% dplyr::filter(get(group) %in% unlist(comparisons)), 
                aes(x=del_length,y=after_stat(density*width))) + 
            geom_histogram(color="black", binwidth=0.2) + 
            scale_x_log10() +
            theme_bw() + xlab("Deletion length (bp)") + ylab("Proportion of deletions") +
            facet_wrap( ~ get(group),scales="fixed")
        ggsave(paste0(outDir,"/deletion_len.histgram_biallelic_vs_control.normalizedInGroup.pdf"),width=10.5,height=5)

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


# compare deletion length for deletions with/without homeology
compare_deletion_length_for_deletions_with_and_without_homeology <- function(deletions, candidates, outDir,
        genotypes = c('BRCA2','control'),
        min_similarity=80, alignment_length_cutoff=30,
        comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control'))
    ){

    library(ggpubr)
    
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    min_similarity_tmp = min_similarity
    alignment_length_cutoff_tmp = alignment_length_cutoff

    dels_homeology <- candidates %>% 
        dplyr::filter(min_similarity==min_similarity_tmp, alignment_length_cutoff==alignment_length_cutoff_tmp) %>% 
        dplyr::mutate(type="With homeology") %>% 
        dplyr::select(SAMPLE.TUMOR, CHROM, start_position, end_position, allelic_status_brca1_brca2, type)

    full_deletions <- deletions %>% 
        left_join(
            get_full_tumor_list() %>%
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>%
            dplyr::select(SAMPLE.TUMOR, allelic_status_brca1_brca2)
        )
    full_deletions <- full_deletions %>% 
        dplyr::left_join(dels_homeology %>% dplyr::select(-allelic_status_brca1_brca2), by=c("SAMPLE.TUMOR","CHROM","start_position","end_position")) %>% 
        dplyr::mutate(type=ifelse(is.na(type),"Without homeology",type)) %>%
        dplyr::filter(allelic_status_brca1_brca2 %in% genotypes) %>%
        dplyr::mutate(group=paste0(type, ' ', allelic_status_brca1_brca2)) %>%
        dplyr::mutate(group=factor(group,levels=unique(unlist(comparisons))))

    fwrite(full_deletions, paste0(outDir,"/full_deletions_with_homeology.tsv"), sep="\t")

    # boxplot deletion length in full deletions and deletions with homeology
    p <- ggplot(full_deletions, 
            aes(x=group,y=del_length)) + 
        geom_violin() +
        geom_boxplot(outlier.shape = NA, width=0.1) + 
        geom_jitter(width=.25, size=0.3, alpha=0.1) + 
        scale_y_log10() +
        stat_compare_means(method="wilcox.test",method.args=list('exact'=FALSE),comparisons=comparisons,label="..p.format..") +
        theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) +
        xlab("") + ylab("Deletion length (bp)")
    ggsave(paste0(outDir,"/deletion_len.boxplot_full_vs_with_homeology.similarity80_hlen30.pdf"),width=5,height=5)

    # histogram deletion length in full deletions and deletions with homeology
    p <- ggplot(full_deletions, 
            aes(x=del_length)) + 
        geom_histogram(color="black", binwidth=0.2) + 
        scale_x_log10() +
        theme_bw() + xlab("Deletion length (bp)") + 
        facet_wrap( ~ group,scales="free_y")
    ggsave(paste0(outDir,"/deletion_len.histogram_full_vs_with_homeology.similarity80_hlen30.pdf"),width=7,height=5)

    

}

make_homeology_rate_plot_in_deletion_sizes <- function(dat, outDir, 
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes.pdf",
    label.y=c(0.3,0.38,0.46), ylim=c(0,0.5)){

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    p <- dat %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            coord_cartesian(ylim = ylim) +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE),
                label.y = label.y, tip.length = 0.01
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
            facet_wrap(~deletion_size,nrow=1) 
    ggsave(paste0(outDir, "/",outfile), width = 7, height = 5)

}


make_homeology_rate_violinplot_in_deletion_sizes <- function(dat, outDir, 
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes.pdf",
    label.y=c(0.3,0.38,0.46), ylim=c(0,0.5)){

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    p <- dat %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = get(group), y = homeology_rate)) +
            geom_violin() +
            geom_sina(size = 0.3) +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            coord_cartesian(ylim = ylim) +
            stat_compare_means(
                comparisons = comparisons,
                method = "wilcox.test", method.args = list("exact" = FALSE),
                label.y = label.y, tip.length = 0.01
            ) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1)) +
            facet_wrap(~deletion_size,nrow=1) 
    ggsave(paste0(outDir, "/",outfile), width = 7, height = 5)

}


make_homeology_rate_plot_in_deletion_sizes_on_x_axis <- function(dat, outDir, 
    group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="make_homeology_rate_plot_in_deletion_sizes_on_x_axis.pdf",
    ylim=c(0,0.5)){

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

    p <- dat %>%
            dplyr::filter(deletions != 0) %>%
            dplyr::filter(get(group) %in% unlist(comparisons)) %>%
            ggplot(aes(x = deletion_size, y = homeology_rate)) +
            geom_boxplot() +
            theme_bw() +
            xlab("") +
            ylab("Homeology rate") +
            coord_cartesian(ylim = ylim) +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
    ggsave(paste0(outDir, "/",outfile), width = 7, height = 5)

}