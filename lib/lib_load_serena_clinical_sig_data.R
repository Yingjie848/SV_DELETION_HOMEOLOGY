# functions to load Serena/Davies' clinical data and siganture data

library(data.table)
library(magrittr)
# library(ipfun)
library(dplyr)


get_serena_signatures <- function() {
    ## Signature and other info from Serena's paper
    serena_sigs <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table21_signature_contribution.csv")
    serena_lst <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table17_HRD_intermediary.csv")
    serena_sigs[, total := rowSums(.SD), .SDcols = 2:13]
    serena_sigs$sig3 <- serena_sigs[[4]] / serena_sigs$total
    # -- For these signature cutoff, for now just using sig3 >=0.3, LST>=20
    serena_sigs$sample <- gsub("(a|b).*$", "", serena_sigs[[1]])
    serena_lst$sample <- gsub("(a|b).*$", "", serena_lst$Sample)
    signatures <- merge(serena_sigs, serena_lst[, .(sample, LST)], by = "sample")
    return(signatures)
}

get_davies_biallelic_samples <- function() {
    sample <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)
    # This is the cohort we are using for final testing 320 ER+/HER2-
    allSamples <- unique(sample$sampleid)

    # Return the biallelic/control samples list based on Davies HRDetect paper
    ## Loading mutation information to determine biallelic vs. control samples
    brca_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4a_mutation_details.csv")
    brca_mut$sample <- gsub("(a|b).*$", "", brca_mut$Sample)
    # -- use this to define the 77 cases with biallelic
    ddr_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4b_otherHRgene_mutation_details.csv")
    ddr_mut$sample <- gsub("(a|b).*$", "", ddr_mut$Sample)

    ddr_mut <- ddr_mut %>% dplyr::filter(Gene %in% c("PALB2","RAD50",'BRIP1')) # -- use these 3 core HRD genes

    # -- use this to remote the other DDR gene (PALB2) mutations
    davies_sample <- fread("input/Davies_HRDetect_2017/Davies_suppltable_1_sampleInfo.csv")

    # Organize breast and OV samples list
    biallelic <- brca_mut[isBrcaMonoallelic == F]
    biallelic_brca1 <- brca_mut[Gene == "BRCA1" & isBrcaMonoallelic == F]
    biallelic_brca2 <- brca_mut[Gene == "BRCA2" & isBrcaMonoallelic == F]
    monoallelic <- ddr_mut[Gene %in% c("PALB2","RAD50","BRIP1")] # -- use these 3 core HRD genes
    signatures <- get_serena_signatures()
    dom_signatures <- signatures[signatures$LST >= 15 & signatures$sig3 >= 0.2]

    # Finalize sample categories
    control_samples <- allSamples[!(allSamples %in% biallelic$sample) &
        !(allSamples %in% monoallelic$sample) &
        !(allSamples %in% dom_signatures$sample)]

    biallelic_samples <- unique(biallelic$sample)
    biallelic_brca1_samples <- unique(biallelic_brca1$sample)
    biallelic_brca2_samples <- unique(biallelic_brca2$sample)
    ddr_samples <- unique(ddr_mut$sample)

    fr <- list()
    fr$control <- control_samples
    fr$biallelic <- biallelic_samples
    fr$biallelic_brca1 <- biallelic_brca1_samples
    fr$biallelic_brca2 <- biallelic_brca2_samples
    fr$ddr_samples <- ddr_samples
    return(fr)
}



get_serena_samples <- function() {
    # Loading Serena breast tumor type info
    sample <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)

    # This is the cohort we are using for final testing 320 ER+/HER2-
    er_pos_sample <- sample[final.ER == "positive" & final.HER2 == "negative"]
    erPosSamples <- unique(er_pos_sample$sampleid)
    trp_neg_sample <- sample[final.ER == "negative" &
        final.HER2 == "negative" &
        final.PR == "negative"]
    trpNegSamples <- unique(trp_neg_sample$sampleid)

    rl <- list()
    rl$erPosSamples <- erPosSamples
    rl$trpNegSamples <- trpNegSamples
    return(rl)
}

add_tumor_type_and_allelic_status <- function(tbl) {


    # Loading existing biallelic information
    # source("lib/lib_load_serena_clinical_sig_data.R")
    sigs <- get_serena_signatures()
    serena <- get_serena_samples()
    allelic_info <- get_davies_biallelic_samples()

    # Use the clinical info from Davies' sample for ER status
    # davies_sample <- fread("input/Davies_HRDetect_2017/Davies_suppltable_1_sampleInfo.csv")
    # names(davies_sample) <- gsub(" ","",colnames(davies_sample))
    # davies_sample$sample <- gsub("(a|b).*$", "", davies_sample$Sample)
    # table(davies_sample$ERstatus)
    # table(davies_sample$isUsedForEvaluation)
    # table(davies_sample$ERstatus, davies_sample$isUsedForEvaluation)

    # Get ER positive samples
    # er <- davies_sample[ERstatus == "positive", ]
    # triple <- davies_sample[ERstatus == "negative", ]
    # -- For ER, triple negative status, using the information from Serena's table (not Davies)

    # Lable ER status and biallelic status
    tbl$tumor_subtype <- ifelse(tbl$SAMPLE.TUMOR %in% serena$erPosSamples, "ER+", "Triple-")
    tbl$allelic_status <- ifelse(tbl$SAMPLE.TUMOR %in% allelic_info$biallelic, "biallelic", ifelse(tbl$SAMPLE.TUMOR %in% allelic_info$control, "control", "other"))
    tbl$allelic_status_brca1_brca2 <- ifelse(tbl$SAMPLE.TUMOR %in% allelic_info$biallelic_brca1, "BRCA1", ifelse(tbl$SAMPLE.TUMOR %in% allelic_info$biallelic_brca2, "BRCA2", ifelse(tbl$SAMPLE.TUMOR %in% allelic_info$control, "control", "other")))

    tbl
}

get_full_tumor_list <- function() {

    # Loading Serena breast tumor type info
    serena <- get_serena_samples()

    # loading sample clinical data
    sample <- fread("input/serena_nature_2016_breast_tumor_560/serena_nature_2016_560_breast_tumor_table1_clindata.csv")
    sample$sampleid <- gsub("(a|b).*$", "", sample$sample_name)

    # load Davies
    allelic_info <- get_davies_biallelic_samples()

    # Lable ER status and biallelic status
    sample$tumor_subtype <- ifelse(sample$sampleid %in% serena$erPosSamples, "ER+", "Triple-")
    sample$allelic_status <- ifelse(sample$sampleid %in% allelic_info$biallelic, "biallelic", ifelse(sample$sampleid %in% allelic_info$control, "control", "other"))
    sample$allelic_status_brca1_brca2 <- ifelse(sample$sampleid %in% allelic_info$biallelic_brca1, "BRCA1",
        ifelse(sample$sampleid %in% allelic_info$biallelic_brca2, "BRCA2",
            ifelse(sample$sampleid %in% allelic_info$control, "control",
                "other"
            )
        )
    )
    sample$other_hrd_genes <- ifelse(sample$sampleid %in% allelic_info$ddr_samples, "Yes", "No")

    # add hrd gene name
    ddr_mut <- fread("input/Davies_HRDetect_2017/Davies_suppltable_4b_otherHRgene_mutation_details.csv")
    ddr_mut$sample <- gsub("(a|b).*$", "", ddr_mut$Sample)

    sample <- sample %>%
        left_join(ddr_mut %>% dplyr::select(sample, Gene) %>% group_by(sample) %>% summarise(Gene=paste0(Gene,collapse=',')),by = c("sampleid"="sample"))

    sample
}

get_hrd <- function(full_tumor_list) {
    hrdetect <- fread("input/Davies_HRDetect_2017/David_2017_HRDetect_score.csv")

    full_tumor_list <- full_tumor_list %>%
        dplyr::mutate(
            other_hrd_genes = ifelse(other_hrd_genes == "Yes", "HRD",
                ifelse(other_hrd_genes == "No", "Non-HRD", "other")
            ),
            hrd_mut = ifelse(allelic_status_brca1_brca2 == "BRCA1", "BRCA1",
                ifelse(allelic_status_brca1_brca2 == "BRCA2", "BRCA2",
                    ifelse(other_hrd_genes == "HRD", "Other HRD",
                        ifelse(allelic_status_brca1_brca2 == "control", "Non-HRD", "other")
                    )
                )
            ),
            hrd_status = ifelse(allelic_status_brca1_brca2 == "BRCA1" | allelic_status_brca1_brca2 == "BRCA2" | other_hrd_genes == "HRD", "HRD",
                ifelse(allelic_status_brca1_brca2 == "control", "Non-HRD", "other")
            )
        ) %>%
        dplyr::mutate(
            hrd_mut = factor(hrd_mut, levels = c("BRCA1", "BRCA2", "Other HRD", "Non-HRD", "other")),
            hrd_status = factor(hrd_status, levels = c("HRD", "Non-HRD", "other"))
        )

    full_tumor_list <- full_tumor_list %>%
        left_join(hrdetect %>%
            dplyr::rename("sample_name" = Sample, "HRDetect_score" = predictorProb) %>%
            dplyr::select(sample_name, HRDetect_score) %>%
            dplyr::mutate(sample_name = gsub("(a|b).*$", "", sample_name)),
        by = "sample_name"
        ) %>%
        dplyr::mutate(
            HRDetect_group = ifelse(HRDetect_score < 0.2, "HRDetect-low",
                ifelse(HRDetect_score >= 0.7, "HRDetect-high", "HRDectect-intermediate")
            )
        ) %>% 
        dplyr::mutate(
            hrd_status_hrdetect = ifelse(hrd_status=='HRD' | HRDetect_group=='HRDetect-high','HRD',
                ifelse(hrd_status=='Non-HRD','Non-HRD','other')
            )
        )
}