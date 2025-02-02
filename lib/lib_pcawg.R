

process_pcawg_sv_deletions <- function(sample_info){

    load_files <- function(dir){
        files <- paste0(dir,"/",list.files(dir,"bedpe.gz$"))
        content <- do.call(rbind,lapply(files,function(x){
            aliquot_id = str_remove(x,".*/")
            aliquot_id = str_remove(aliquot_id,".pcawg.*$")
            out = fread(x)
            out$tumour_specimen_aliquot_id = aliquot_id
            out
        }))
    }

    if(!file.exists("data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv.icgc.rds") | !file.exists("data/PCAWG_consensus_sv/PCAWG_consensus_sv.tcga.rds")){
        consensus_sv_icgc <- load_files(dir="data/pcawg/PCAWG_consensus_sv/icgc/open/")
        consensus_sv_tcga <- load_files(dir="data/pcawg/PCAWG_consensus_sv/tcga/open/")

        length(unique(consensus_sv_icgc$tumour_specimen_aliquot_id))
        length(unique(consensus_sv_tcga$tumour_specimen_aliquot_id))

        saveRDS(consensus_sv_icgc,file="data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv.icgc.rds")
        saveRDS(consensus_sv_tcga,file="data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv.tcga.rds")
    }else{
        print("Loading existing PCAWG_consensus_sv.icgc.rds and PCAWG_consensus_sv.tcga.rds")
        consensus_sv_icgc = readRDS("data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv.icgc.rds")
        consensus_sv_tcga = readRDS("data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv.tcga.rds")
    }

    # count number of deletions
    sv_dels_icgc <- consensus_sv_icgc %>% 
        dplyr::filter(svclass=='DEL') %>%
        dplyr::mutate(del_length=start2-end1)
    sv_dels_icgc %>% 
        fwrite("data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv_deletions.icgc.tsv",sep="\t")

    sv_dels_tcga <- consensus_sv_tcga %>% 
        dplyr::filter(svclass=='DEL') %>%
        dplyr::mutate(del_length=start2-end1)
    sv_dels_tcga %>% 
        fwrite("data/pcawg/PCAWG_consensus_sv/PCAWG_consensus_sv_deletions.tcga.tsv",sep="\t")

}



summarize_candidates_by_dcc_project_code <- function(candidates, full_deletions, outDir, min_number_of_deletions_per_sample=0) {
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
            n_events = length(SAMPLE.TUMOR),
            mean_del_length = mean(del_length),
            mean_alignment_length = mean(alignment_length),
            mean_similarity = mean(identity)
        )
    sample_summary %>% fwrite(paste0(outDir, "/sample_summary.tsv"), sep = "\t")

    # get full tumor sample list, add no. of events per sample, total number of deletions per sample
    deletions_per_sample <- full_deletions %>%
        group_by(SAMPLE.TUMOR,dcc_project_code,project_code) %>%
        dplyr::summarise(deletions = length(SAMPLE.TUMOR))
    sample <- deletions_per_sample %>%
        left_join(sample_summary, by = "SAMPLE.TUMOR") %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events),
            deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        ) %>%
        dplyr::filter(deletions >= min_number_of_deletions_per_sample)
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

    # summarize samples by project_code, calculate proportion of homeology samples, proportion of homeology events
    homeology_rate_per_group <- sample %>%
        group_by(dcc_project_code) %>%
        summarise(total_samples = length(unique(SAMPLE.TUMOR)), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE), median_homeology_rate=median(homeology_rate)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        dplyr::arrange(desc(median_homeology_rate))
        
    homeology_rate_per_group %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumor_types.tsv"), sep = "\t")

    # make homeology_rate (n_events/deletions) boxplot/violin by dcc_project_code
    p <- sample %>% 
        dplyr::mutate(dcc_project_code=factor(dcc_project_code,levels=homeology_rate_per_group$dcc_project_code)) %>%
        dplyr::filter(deletions != 0, !is.na(dcc_project_code), dcc_project_code %in% homeology_rate_per_group$dcc_project_code) %>%
        ggplot(aes(x = dcc_project_code, y = homeology_rate,color=project_code)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none")
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.pdf"), width = 10, height = 5)

    p <- sample %>% 
        dplyr::mutate(dcc_project_code=factor(dcc_project_code,levels=homeology_rate_per_group$dcc_project_code)) %>%
        dplyr::filter(deletions != 0, !is.na(dcc_project_code)) %>%
        ggplot(aes(x = dcc_project_code, y = homeology_rate,color=project_code)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none") +
        coord_cartesian(ylim=c(0,0.2))
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.zoom.pdf"), width = 10, height = 5)


    # for samples with deletions

    sample <- sample %>% dplyr::filter(deletions > 0, deletions >= min_number_of_deletions_per_sample)

    sample %>%
        summarise(
            total_samples = length(SAMPLE.TUMOR),
            total_deletions = sum(deletions),
            homeology_samples = sum(n_events >= 1, na.rm = TRUE),
            total_homeology_events = sum(n_events, na.rm = TRUE)
        ) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumors_with_deletions.tsv"), sep = "\t")

    # make plots
    make_deletion_level_plots(candidates = candidates, outDir = outDir)

}


summarize_candidates_by_project_code <- function(candidates, full_deletions, outDir, min_number_of_deletions_per_sample=0) {

    dir.create(outDir)
    
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
            n_events = length(SAMPLE.TUMOR),
            mean_del_length = mean(del_length),
            mean_alignment_length = mean(alignment_length),
            mean_similarity = mean(identity)
        )
    sample_summary %>% fwrite(paste0(outDir, "/sample_summary.tsv"), sep = "\t")

    # get full tumor sample list, add no. of events per sample, total number of deletions per sample
    deletions_per_sample <- full_deletions %>%
        group_by(SAMPLE.TUMOR,project_code) %>%
        dplyr::summarise(deletions = length(SAMPLE.TUMOR))
    sample <- deletions_per_sample %>%
        left_join(sample_summary, by = "SAMPLE.TUMOR") %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events),
            deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        ) %>%
        dplyr::filter(deletions >= min_number_of_deletions_per_sample)
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

    # summarize samples by project_code, calculate proportion of homeology samples, proportion of homeology events
    homeology_rate_per_group <- sample %>%
        group_by(project_code) %>%
        summarise(total_samples = length(unique(SAMPLE.TUMOR)), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE), median_homeology_rate=median(homeology_rate)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        dplyr::arrange(desc(median_homeology_rate))
        
    homeology_rate_per_group %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumor_types.tsv"), sep = "\t")

    # exclude project_code with <10 samples
    projects_with_enough_samples <- homeology_rate_per_group %>% dplyr::filter(total_samples >= 10) %>% dplyr::pull(project_code)
    sample <- sample %>% dplyr::filter(project_code %in% projects_with_enough_samples)

    # make homeology_rate (n_events/deletions) boxplot/violin by project_code
    p <- sample %>% 
        dplyr::mutate(project_code=factor(project_code,levels=homeology_rate_per_group$project_code)) %>%
        dplyr::filter(deletions != 0, !is.na(project_code), project_code %in% homeology_rate_per_group$project_code) %>%
        ggplot(aes(x = project_code, y = homeology_rate,color=project_code)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none")
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.pdf"), width = 10, height = 5)

    p <- sample %>% 
        dplyr::mutate(project_code=factor(project_code,levels=homeology_rate_per_group$project_code)) %>%
        dplyr::filter(deletions != 0, !is.na(project_code)) %>%
        ggplot(aes(x = project_code, y = homeology_rate,color=project_code)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none") +
        coord_cartesian(ylim=c(0,0.2))
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.zoom.pdf"), width = 10, height = 5)


    # for samples with deletions

    sample <- sample %>% dplyr::filter(deletions > 0, deletions >= min_number_of_deletions_per_sample)

    sample %>%
        summarise(
            total_samples = length(SAMPLE.TUMOR),
            total_deletions = sum(deletions),
            homeology_samples = sum(n_events >= 1, na.rm = TRUE),
            total_homeology_events = sum(n_events, na.rm = TRUE)
        ) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumors_with_deletions.tsv"), sep = "\t")

    # make plots
    make_deletion_level_plots(candidates = candidates, outDir = outDir)

}


summarize_candidates_by_dcc_project_code_hpv <- function(candidates, full_deletions, outDir, min_number_of_deletions_per_sample=0) {

    dir.create(outDir)
    
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
            n_events = length(SAMPLE.TUMOR),
            mean_del_length = mean(del_length),
            mean_alignment_length = mean(alignment_length),
            mean_similarity = mean(identity)
        )
    sample_summary %>% fwrite(paste0(outDir, "/sample_summary.tsv"), sep = "\t")

    # get full tumor sample list, add no. of events per sample, total number of deletions per sample
    deletions_per_sample <- full_deletions %>%
        group_by(SAMPLE.TUMOR,dcc_project_code_hpv) %>%
        dplyr::summarise(deletions = length(SAMPLE.TUMOR))
    sample <- deletions_per_sample %>%
        left_join(sample_summary, by = "SAMPLE.TUMOR") %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events),
            deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        ) %>%
        dplyr::filter(deletions >= min_number_of_deletions_per_sample)
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

    # summarize samples by dcc_project_code_hpv, calculate proportion of homeology samples, proportion of homeology events
    homeology_rate_per_group <- sample %>%
        group_by(dcc_project_code_hpv) %>%
        summarise(total_samples = length(unique(SAMPLE.TUMOR)), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE), median_homeology_rate=median(homeology_rate)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        dplyr::arrange(desc(median_homeology_rate))
        
    homeology_rate_per_group %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumor_types.tsv"), sep = "\t")

    # make homeology_rate (n_events/deletions) boxplot/violin by dcc_project_code_hpv
    p <- sample %>% 
        dplyr::mutate(dcc_project_code_hpv=factor(dcc_project_code_hpv,levels=homeology_rate_per_group$dcc_project_code_hpv)) %>%
        dplyr::filter(deletions != 0, !is.na(dcc_project_code_hpv), dcc_project_code_hpv %in% homeology_rate_per_group$dcc_project_code_hpv) %>%
        ggplot(aes(x = dcc_project_code_hpv, y = homeology_rate,color=dcc_project_code_hpv)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none")
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.pdf"), width = 10, height = 5)

    p <- sample %>% 
        dplyr::mutate(dcc_project_code_hpv=factor(dcc_project_code_hpv,levels=homeology_rate_per_group$dcc_project_code_hpv)) %>%
        dplyr::filter(deletions != 0, !is.na(dcc_project_code_hpv)) %>%
        ggplot(aes(x = dcc_project_code_hpv, y = homeology_rate,color=dcc_project_code_hpv)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none") +
        coord_cartesian(ylim=c(0,0.2))
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.homeology_rate.zoom.pdf"), width = 10, height = 5)


    # make deletion boxplot/violin by dcc_project_code_hpv
    p <- sample %>% 
        dplyr::mutate(dcc_project_code_hpv=factor(dcc_project_code_hpv,levels=homeology_rate_per_group$dcc_project_code_hpv)) %>%
        dplyr::filter(!is.na(dcc_project_code_hpv), dcc_project_code_hpv %in% homeology_rate_per_group$dcc_project_code_hpv) %>%
        ggplot(aes(x = dcc_project_code_hpv, y = deletions,color=dcc_project_code_hpv)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Deletions") +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none")
    ggsave(paste0(outDir, "/ssa_events_in_tumor_types.deletions.pdf"), width = 10, height = 5)


    # make homeology_rate (n_events/deletions) boxplot/violin for HNSC-US_HPV+, HNSC-US_HPV-, CESC-US
    p <- sample %>% 
        dplyr::filter(dcc_project_code_hpv %in% c('HNSC-US_HPV+','HNSC-US_HPV-','CESC-US')) %>% 
        dplyr::mutate(dcc_project_code_hpv=factor(dcc_project_code_hpv,levels=c('HNSC-US_HPV+','HNSC-US_HPV-','CESC-US'))) %>%
        dplyr::filter(deletions != 0) %>%
        ggplot(aes(x = dcc_project_code_hpv, y = homeology_rate,color=dcc_project_code_hpv)) +
        geom_boxplot(outlier.shape=NA) +
        theme_bw() +
        geom_jitter(width=0.2,size=0.1) +
        xlab("") +
        ylab("Homeology rate") +
        stat_compare_means(comparisons=list(c('HNSC-US_HPV+','HNSC-US_HPV-'),c('HNSC-US_HPV+','CESC-US'),c('HNSC-US_HPV-','CESC-US'))) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1),legend.position = "none")
    ggsave(paste0(outDir, "/ssa_events_in_HNSC_CESC.homeology_rate.pdf"), width = 5, height = 5)


    # for samples with deletions

    sample <- sample %>% dplyr::filter(deletions > 0, deletions >= min_number_of_deletions_per_sample)

    sample %>%
        summarise(
            total_samples = length(SAMPLE.TUMOR),
            total_deletions = sum(deletions),
            homeology_samples = sum(n_events >= 1, na.rm = TRUE),
            total_homeology_events = sum(n_events, na.rm = TRUE)
        ) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_tumors_with_deletions.tsv"), sep = "\t")

    # make plots
    make_deletion_level_plots(candidates = candidates, outDir = outDir)

}


get_candidates <- function(full_deletions, blast_output_same_strand, outDir, combineOnly = FALSE, take_longest = TRUE, take_most_similar = FALSE, min_number_of_deletions_per_sample = 0) {

    #
    similarities <- c(0, 80, 90, 100)
    alignment_lens <- c(5, 10, 15, 20, 25, 30, 35, 40)

    if (!combineOnly) {
        for (min_similarity in similarities) {
            cat("Similarity:", min_similarity, "\n")
            for (alignment_length_cutoff in alignment_lens) {
                cat("  Min alignment length:", alignment_length_cutoff, "\n")

                outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")
                dir.create(outdir_candidates)

                # for each event, get longest alignment with >=X% similarity
                candidates <- identify_ssa_candidates(blast_output_same_strand, min_similarity = min_similarity, alignment_length_cutoff = alignment_length_cutoff, take_longest = take_longest, take_most_similar = take_most_similar)
                # add deletion length
                candidates <- candidates %>% left_join(full_deletions, by = c("SAMPLE.TUMOR", "CHROM", "start_position", "end_position"))

                candidates %>% fwrite(paste0(outdir_candidates, "/candidates.tsv"), sep = "\t")

                if (nrow(candidates) > 0) {
                    summarize_candidates_by_dcc_project_code(candidates = candidates, full_deletions = full_deletions, outDir = outdir_candidates, min_number_of_deletions_per_sample = min_number_of_deletions_per_sample)
                }
            }
        }


        # combine candidates.tsv and ssa_events_for_each_sample.tsv for different parameters
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
    }

    list(candidates, sample)
}



make_boxplot_by_project_code <- function(full_deletions,outDir, min_number_of_deletions_per_sample = 0) {

    #
    similarities <- c(0, 80, 90, 100)
    alignment_lens <- c(5, 10, 15, 20, 25, 30, 35, 40)

    for (min_similarity in similarities) {
        cat("Similarity:", min_similarity, "\n")
        for (alignment_length_cutoff in alignment_lens) {
            cat("  Min alignment length:", alignment_length_cutoff, "\n")

            outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")
            dir.create(outdir_candidates)

            candidates <- fread(paste0(outdir_candidates, "/candidates.tsv"))

            if (nrow(candidates) > 0) {
                summarize_candidates_by_project_code(candidates = candidates, full_deletions = full_deletions, outDir = paste0(outdir_candidates,"/plots_by_project_code"), min_number_of_deletions_per_sample = min_number_of_deletions_per_sample)
            }
        }
    }
}

make_boxplot_for_dcc_project_code_hpv <- function(full_deletions,outDir, min_number_of_deletions_per_sample = 0, 
    similarities = c(0, 80, 90, 100), 
    alignment_lens = c(5, 10, 15, 20, 25, 30, 35, 40)) {

    for (min_similarity in similarities) {
        cat("Similarity:", min_similarity, "\n")
        for (alignment_length_cutoff in alignment_lens) {
            cat("  Min alignment length:", alignment_length_cutoff, "\n")

            outdir_candidates <- paste0(outDir, "/results_similarity_", min_similarity, "pct_", alignment_length_cutoff, "bp")
            dir.create(outdir_candidates)

            candidates <- fread(paste0(outdir_candidates, "/candidates.tsv"))

            if (nrow(candidates) > 0) {
                summarize_candidates_by_dcc_project_code_hpv(candidates = candidates, full_deletions = full_deletions, outDir = paste0(outdir_candidates,"/plots_by_dcc_project_code_hpv"), min_number_of_deletions_per_sample = min_number_of_deletions_per_sample)
            }
        }
    }
}