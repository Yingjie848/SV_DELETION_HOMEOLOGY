# take blast table and identify hits most likely related to SSA

# take longest alignment with at least X% similarity
# or take shortest alignment with at least X% similarity if take_longest=FALSE
identify_ssa_candidates <- function(bltbl, min_similarity = 0, alignment_length_cutoff = 0, take_longest = TRUE, take_most_similar = FALSE) {

    bltbl <- bltbl %>% dplyr::filter(q_strand == "+", s_strand == "+")

    if (take_longest) { # take the longest alignment
        candidates <- bltbl %>%
            dplyr::filter(identity >= min_similarity, alignment_length >= alignment_length_cutoff) %>%
            dplyr::group_by(SAMPLE.TUMOR, CHROM, start_position, end_position) %>%
            dplyr::arrange(desc(alignment_length)) %>%
            dplyr::slice_head(n = 1)
    } else {
        if (take_most_similar) {
            # take the most similar
            candidates <- bltbl %>%
                dplyr::filter(identity >= min_similarity, alignment_length >= alignment_length_cutoff) %>%
                dplyr::group_by(SAMPLE.TUMOR, CHROM, start_position, end_position) %>%
                dplyr::arrange(desc(identity), alignment_length) %>%
                dplyr::slice_head(n = 1)
        } else {
            # take the shortest alignment
            candidates <- bltbl %>%
                dplyr::filter(identity >= min_similarity, alignment_length >= alignment_length_cutoff) %>%
                dplyr::group_by(SAMPLE.TUMOR, CHROM, start_position, end_position) %>%
                dplyr::arrange(alignment_length) %>%
                dplyr::slice_head(n = 1)
        }
    }
}

# summarize candidates by sample
summarize_candidates <- function(candidates, d.out, outDir) {
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
            tumor_subtype = tumor_subtype[1],
            allelic_status = allelic_status[1],
            allelic_status_brca1_brca2 = allelic_status_brca1_brca2[1],
            n_events = length(SAMPLE.TUMOR),
            mean_del_length = mean(del_length),
            mean_alignment_length = mean(alignment_length),
            mean_similarity = mean(identity)
        )
    sample_summary %>% fwrite(paste0(outDir, "/sample_summary.tsv"), sep = "\t")

    # get full tumor sample list, add no of events per sample, total number of deletions per sample
    deletions_per_sample <- d.out %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(deletions = length(SAMPLE.TUMOR))
    sample <- get_full_tumor_list() %>%
        get_hrd() %>%
        merge(., sample_summary %>% dplyr::select(-tumor_subtype, -allelic_status, -allelic_status_brca1_brca2), by.x = "sampleid", by.y = "SAMPLE.TUMOR", all.x = TRUE) %>%
        merge(., deletions_per_sample, by.x = "sampleid", by.y = "SAMPLE.TUMOR", all.x = TRUE) %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events), deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        )
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

    # summarize samples by BRCA1+BRCA2/control, calculate proportion of homeology samples, proportion of homeology events
    sample %>%
        group_by(allelic_status) %>%
        summarise(total_samples = length(allelic_status), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_bialliec_tumors.tsv"), sep = "\t")

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

    # make homeology_rate (n_events/deletions) boxplot by BRCA1/BRCA2/control
    p <- sample %>%
        dplyr::filter(deletions != 0, allelic_status_brca1_brca2 != "other") %>%
        ggplot(aes(x = allelic_status_brca1_brca2, y = homeology_rate)) +
        geom_boxplot() +
        theme_bw() +
        geom_jitter() +
        xlab("") +
        ylab("Homeology rate") +
        stat_compare_means(
            comparisons = list(c("BRCA1", "BRCA2"), c("BRCA1", "control"), c("BRCA2", "control")),
            method = "wilcox.test", method.args = list("exact" = FALSE)
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
    ggsave(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.pdf"), width = 3, height = 5)

    p <- sample %>%
        dplyr::filter(deletions != 0, allelic_status_brca1_brca2 != "other") %>%
        ggplot(aes(x = allelic_status_brca1_brca2, y = homeology_rate)) +
        geom_violin(position = "identity") +
        geom_sina(size = 0.3) +
        theme_bw() +
        xlab("") +
        ylab("Homeology rate") +
        stat_compare_means(
            comparisons = list(c("BRCA1", "BRCA2"), c("BRCA1", "control"), c("BRCA2", "control")),
            method = "wilcox.test", method.args = list("exact" = FALSE)
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))
    ggsave(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors.homeology_rate.violin.pdf"), width = 3, height = 5)


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
    ggsave(paste0(outDir, "/ssa_events_in_hrd_tumors.homeology_rate.violin.pdf"), width = 3, height = 5)

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
    ggsave(paste0(outDir, "/ssa_events_in_other_hrd_tumors.homeology_rate.violin.pdf"), width = 3, height = 5)


    # for samples with deletions

    sample <- sample %>% dplyr::filter(deletions > 0)

    # summarize samples by BRCA1+BRCA2/control, calculate proportion of homeology samples, proportion of homeology events
    sample %>%
        group_by(allelic_status) %>%
        summarise(total_samples = length(allelic_status), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_bialliec_tumors_with_deletions.tsv"), sep = "\t")

    # summarize samples by BRCA1/BRCA2/control, calculate proportion of homeology samples, proportion of homeology events
    sample %>%
        group_by(allelic_status_brca1_brca2) %>%
        summarise(total_samples = length(allelic_status_brca1_brca2), total_deletions = sum(deletions), homeology_samples = sum(n_events >= 1, na.rm = TRUE), total_homeology_events = sum(n_events, na.rm = TRUE)) %>%
        dplyr::mutate(
            proportion_of_homeology_samples = homeology_samples / total_samples,
            mean_homeology_events_per_sample = total_homeology_events / total_samples,
            proportion_of_homeology_events = total_homeology_events / total_deletions
        ) %>%
        fwrite(paste0(outDir, "/ssa_events_in_bialliec_brca1_brca2_tumors_with_deletions.tsv"), sep = "\t")

    # make plots
    analyze_ssa_homeology_events(candidates = candidates, outDir = outDir)

    # analyze_ssa_homeology_events_by_tumor_subtype_and_allelic_status(candidates=candidates,outDir=outDir)

    # analyze_ssa_homeology_events_by_tumor_subtype_and_allelic_status_averaged_by_tumor(candidates=candidates,outDir=outDir)
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
