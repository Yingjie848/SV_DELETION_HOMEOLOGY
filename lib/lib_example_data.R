
# summarize candidates by sample
summarize_candidates <- function(candidates, full_deletions, outDir) {

    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)


    # make sample level summary

    # get full tumor sample list, add no of events per sample, total number of deletions per sample, etc.
    sample_summary <- candidates %>%
        group_by(SAMPLE.TUMOR) %>%
        dplyr::summarise(
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
    sample <- 
        merge(sample_summary, deletions_per_sample, by.x = "SAMPLE.TUMOR", by.y = "SAMPLE.TUMOR", all.x = TRUE) %>%
        dplyr::mutate(
            n_events = ifelse(is.na(n_events), 0, n_events), 
            deletions = ifelse(is.na(deletions), 0, deletions),
            homeology_rate = n_events / deletions
        )
    sample %>% fwrite(paste0(outDir, "/ssa_events_for_each_sample.tsv"), sep = "\t")

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
                #candidates <- add_tumor_type_and_allelic_status(candidates)
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

