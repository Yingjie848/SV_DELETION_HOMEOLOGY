rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(tidyr)

setwd("/home/zhuy1/my_projects_nrlab/Manisha_SSA/ssa_homeology_20231119")

source("lib/lib_load_serena_clinical_sig_data.R")
source("lib/lib_blast.R")
source("lib/lib_identify_homeology.R")
source("lib/lib_plots.R")
source("lib/lib_repeats.R")
source("lib/lib_serena_breast_data.R")


outdir="output/serena_SV_deletions_ssa_blast_lr100bp"; dir.create(outdir)
blast_dir=paste0(outdir,"/blast_output"); dir.create(blast_dir)

hrdetect <- fread("input/Davies_HRDetect_2017/David_2017_HRDetect_score.csv")


## -----------------------------------------------------------------------------------------------
# prepare SV deletions from Serena data

prepare_data <- function(){

    # Read in the data
    d <- fread("data/serena_data/BRCA-EU_SV/structural_somatic_mutation.BRCA-EU.tsv")
    table(d$variant_type)

    # Filter for deletions
    d <- d[variant_type == "deletion"]

    ## Snapshot of the deletion data
    # Classes ‘data.table’ and 'data.frame':  17564 obs. of  47 variables:
    #  $ icgc_donor_id                   : chr  "DO218489" "DO218489" "DO218489" "DO218489" ...
    #  $ project_code                    : chr  "BRCA-EU" "BRCA-EU" "BRCA-EU" "BRCA-EU" ...
    #  $ icgc_specimen_id                : chr  "SP117710" "SP117710" "SP117710" "SP117710" ...
    #  $ icgc_sample_id                  : chr  "SA543682" "SA543682" "SA543682" "SA543682" ...
    #  $ submitted_sample_id             : chr  "PD8623a" "PD8623a" "PD8623a" "PD8623a" ...
    #  $ submitted_matched_sample_id     : chr  "PD8623b" "PD8623b" "PD8623b" "PD8623b" ...
    #  $ variant_type                    : chr  "deletion" "deletion" "deletion" "deletion" ...
    #  $ sv_id                           : chr  "CGP_basis_stsm_7558" "CGP_basis_stsm_7545" "CGP_basis_stsm_7541" "CGP_basis_stsm_7551" ...
    #  $ placement                       : int  1 1 1 1 1 1 1 1 1 1 ...
    #  $ annotation                      : chr  "Chr.11  133629242(46)--AGAA--133643401(05)  Chr.11  (score 97)" "Chr.3  13256939][13260353  Chr.3  (score 98)" "Chr.1  75144883(89)--TTCTGC--104591725(31)  Chr.1  (score 97)" "Chr.5  135585320(21)--A--135866679(80)  Chr.5  (score 100)" ...
    #  $ interpreted_annotation          : logi  NA NA NA NA NA NA ...
    #  $ chr_from                        : chr  "11" "3" "1" "5" ...
    #  $ chr_from_bkpt                   : int  133629244 13256939 75144886 135585320 11983732 12457507 3309561 76658015 55568194 108189931 ...
    #  $ chr_from_strand                 : int  1 1 1 1 1 1 1 1 1 1 ...
    #  $ chr_from_range                  : int  4 0 6 1 1 0 2 1 1 2 ...
    #  $ chr_from_flanking_seq           : logi  NA NA NA NA NA NA ...
    #  $ chr_to                          : chr  "11" "3" "1" "5" ...
    #  $ chr_to_bkpt                     : int  133643403 13260353 104591728 135866679 29549859 12459976 3415867 76664753 56161643 109861688 ...
    #  $ chr_to_strand                   : int  1 1 1 1 1 1 1 1 1 1 ...
    #  $ chr_to_range                    : int  4 0 6 1 1 0 2 1 1 2 ...
    #  $ chr_to_flanking_seq             : logi  NA NA NA NA NA NA ...
    #  $ assembly_version                : chr  "GRCh37" "GRCh37" "GRCh37" "GRCh37" ...
    #  $ sequencing_strategy             : chr  "WGS" "WGS" "WGS" "WGS" ...
    #  $ microhomology_sequence          : chr  "AGAA" "" "TTCTGC" "A" ...

    # Encoding fields for the pipeline calculation
    d$sampleid <- gsub("(a|b).*$", "", d$submitted_sample_id) # There are 569 unique sample_id, after removing a|b, it's still 569 unique sample id
    # table(d$verification_status) # all the SV are not tested

    d$SAMPLE.TUMOR  <- d$sampleid
    d$SAMPLE.NORMAL <- d$submitted_matched_sample_id
    d$CHROM         <- d$chr_from


    ##
    # There is a cavet with the exact position of the breakpoints
    # From the input data, it seems that the breakpoints have a range within certian window
    # -- the given microhomology sequence is the exact maching sequence but not exactly at breakpoints
    # -- Looking for the annotation, the start is from breakpoints - from range
    #    the end is to breakpoints + to range
    # d$start_position <- d$chr_from_bkpt - round(d$chr_from_range / 2)
    # d$end_position <- d$chr_to_bkpt - round(d$chr_to_range / 2)
    # -- Current decision: using the given breakpoints (omitting the range info)
    d$start_position <- d$chr_from_bkpt
    d$end_position <- d$chr_to_bkpt
    d$del_length <- d$end_position - d$start_position + 1
    table(d$end_position > d$start_position)
    summary(d$del_length)

    #
    d$Type <- "SOMATIC"
    # Coding variantCaller, if Pindel in the string, then it's Pindel, otherwise it's PCAWG_consensus
    table(d$variation_calling_algorithm)
    d$variantCaller <- "BRASS"

    d.out <- d[, .(
        SAMPLE.TUMOR,
        SAMPLE.NORMAL,
        CHROM,
        start_position,
        end_position,
        del_length,
        variant_type,
        variantCaller,
        assembly_version,
        microhomology_sequence,
        non_templated_sequence
    )] %>% distinct

    # Write out the data
    d.out %>% fwrite(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv"), sep = "\t")

    # make deletion length distribution
    p <- ggplot(d.out,aes(x=del_length)) + 
            geom_histogram(bins=50,fill="gray",color="blue") + 
            scale_x_log10(breaks=c(1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)) + 
            xlab("Deletion length") +
            facet_zoom(xlim = c(100, 50000), zoom.size = 1,split=TRUE)
    ggsave(paste0(outdir,"/serena_BRCA-EU_SV_deletions.deletion_length_histogram.pdf"),width=7,height=7)

    # count number of deletions per sample
    deletions_per_sample <- d.out %>% group_by(SAMPLE.TUMOR) %>% dplyr::summarise(deletions=length(SAMPLE.TUMOR))
    deletions_per_sample %>% fwrite(paste0(outdir,"/serena_BRCA-EU_SV_deletions.deletions_per_sample.tsv"),sep="\t")

    d.out
}

d.out <- prepare_data()


## -----------------------------------------------------------------------------------------------
# Run blastn, aligning 100bp leftside and 100bp rightside of the first breakpoint to  100bp leftside and 100bp rightside of the second breakpoint

# run pairwise blast for two sequences of each deletion
# if everything runs well, delete tmp_dir afterwards, which contains intermediate blast output files
run_pairwise_blast_lr100(
    deletion_tbl=d.out,
    bltbl_file=paste0(blast_dir,"/blast_output.tsv"),
    tmp_dir=paste0(blast_dir,"/tmp"))

# Process blast output, extract alignments with the same direction
# Output: blast_output.same_direction.tsv
process_blast_output_serena_data(blast_dir)


## -----------------------------------------------------------------------------------------------
# identify deletions with homeology candidates by taking most similar and longest alignment
# At the end, we decide to use >=80% similarity, and >=30bp alignment length

##
## for deletions >1kb 
##

outDir=paste0(outdir,"/output_deletions_above_1KB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

## --- extract deletions
deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1000)

## --- count deletions in BRCA2 and control
dels <- deletions %>% 
    left_join(get_full_tumor_list() %>% dplyr::select(sampleid,allelic_status_brca1_brca2),by=c('SAMPLE.TUMOR'='sampleid')) %>%
    dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 11912

# according to Serena's paper, they grouped deletions into 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
dels %>% filter(del_length>1e3,del_length<=1e4) %>% nrow # 2639
dels %>% filter(del_length>1e4,del_length<=1e5) %>% nrow # 1854
dels %>% filter(del_length>1e5,del_length<=1e6) %>% nrow # 2077
dels %>% filter(del_length>1e6,del_length<=1e7) %>% nrow # 2948
dels %>% filter(del_length>1e7) %>% nrow # 2394

## --- count deletions in BRCA1, BRCA2 and control
# according to Serena's paper, they grouped deletions into 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
deletions %>% filter(del_length>1e3,del_length<=1e4) %>% nrow # 4816
deletions %>% filter(del_length>1e4,del_length<=1e5) %>% nrow # 3450
deletions %>% filter(del_length>1e5,del_length<=1e6) %>% nrow # 2918
deletions %>% filter(del_length>1e6,del_length<=1e7) %>% nrow # 3418
deletions %>% filter(del_length>1e7) %>% nrow # 2960

## --- check genotype samples
deletions %>%
    left_join(get_full_tumor_list(),by=c('SAMPLE.TUMOR'='sampleid')) %>%
    group_by(allelic_status_brca1_brca2) %>%
    summarise(samples = length(unique(SAMPLE.TUMOR)))
#  allelic_status_brca1_brca2 samples
#1 BRCA1                           22
#2 BRCA2                           11
#3 control                         62
#4 other                           137

## --- examine deletion length for BRCA1/2/control tumors
examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

## --- examine deletion length for BRCA2/control tumors
examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA2','control')))

## --- load blast output data with the same direction, add deletion length, and filter out deletions without length
blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>%
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length)

## --- get deletions with homeology candidates by similarity + homeology length
#out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))

# For the final analysis, we only use 80% similarity and 30bp alignment length
out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(80), alignment_lens = c(30))

candidates <- out[[1]] # deletions with homeology candidates
sample     <- out[[2]] # sample level results, including homeology rate (proportion of deletions with homeology)

candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")


## --- further analysis

## --- Grid search similarity x homeology seqeunce length
make_sample_plots_similarity_vs_alignment_length(sample=sample,outDir=outDir)

# for 80% similaritys
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80"),
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

# for 80% similarity, BRCA2/control
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'),min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80_BRCA2_control"),
    comparisons=list(c('BRCA2','control')))

## --- check controls with long tails, it could be caused by low number of deletions
p <- sample %>% dplyr::filter(min_similarity==80,alignment_length_cutoff==30,allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>%
    ggplot(aes(x=deletions,y=homeology_rate)) + geom_point() + scale_x_log10() + facet_wrap(~allelic_status_brca1_brca2,nrow=1)
dir.create(paste0(outDir,"/plots_homeology_rate_vs_deletions"))
ggsave(paste0(outDir,"/plots_homeology_rate_vs_deletions/homeology_rate_vs_deletions.pdf"),p,width=10,height=5)


## --- compare alignment length for most similar and longest alignment
compare_alignment_length_most_similar_and_longest(candidates,outDir)


## --- check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
make_distance_to_breakpoint_plot(candidates,outDir=outDir,similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))


## --- analyze repeats

# for candidates
candidates <- fread(paste0(outDir,"/ssa_events_candidates.tsv"))
candidates <- annotate_homeology_sequence_to_repeats(candidates)
candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.add_repeats.tsv"),sep="\t")

summarize_repeats(candidates,outDir=paste0(outDir,"/results_repeats"))
summarize_repeats_for_deletions(candidates %>% filter(min_similarity==80,alignment_length_cutoff==30),outDir=paste0(outDir,"/results_repeats/deletions_with_homeology"))

# for all deletions
deletions_add_repeats <- annotate_breakpoint_flanking_sequence_to_repeats(deletions)
deletions_add_repeats %>% fwrite(paste0(outDir,"/serena_BRCA-EU_SV_deletions.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))


## --- compare deletion length for deletions with/without homeology
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_30bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=30,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))
    
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_25bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=25,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))


## -- check POLQ expression in HRDetect high and low groups
check_POLQ_expression_in_HRDetect_group()


##
## for deletions 1kb-100kb
##

outDir=paste0(outdir,"/output_deletions_1KB-100KB_final"); dir.create(outDir)

## --- extract deletions
deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e3, del_length<1e5)

## --- count deletions in BRCA2 and control
dels <- deletions %>% 
    left_join(get_full_tumor_list() %>% dplyr::select(sampleid,allelic_status_brca1_brca2),by=c('SAMPLE.TUMOR'='sampleid')) %>%
    dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 4493

## --- check genotype samples
deletions %>%
    left_join(get_full_tumor_list(),by=c('SAMPLE.TUMOR'='sampleid')) %>%
    group_by(allelic_status_brca1_brca2) %>%
    summarise(samples = length(unique(SAMPLE.TUMOR)))
#allelic_status_brca1_brca2 samples
#  <chr>                        <int>
#1 BRCA1                           47
#2 BRCA2                           30
#3 control                        335
#4 other                           55

## --- examine deletion length for BRCA1/2/control tumors
examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

## --- examine deletion length for BRCA2/control tumors
examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA2','control')))

# load blast output data with the same direction, add deletion length, and filter out deletions without length
blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>%
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length)

## --- get SSA homeology candidates by similarity + homeology length
#out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))
# for now, we only use 80% similarity and 30bp alignment length
out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(80), alignment_lens = c(30))

candidates <- out[[1]]
sample     <- out[[2]]

candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

## --- further analysis

## --- Grid search similarity x homeology seqeunce length
make_sample_plots_similarity_vs_alignment_length(sample=sample,outDir=outDir)


# for 80% similaritys
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80"),
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

# for 80% similarity, BRCA2/control
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'),min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80_BRCA2_control"),
    comparisons=list(c('BRCA2','control')))

## --- check controls with long tails, it could be caused by low number of deletions
p <- sample %>% dplyr::filter(min_similarity==80,alignment_length_cutoff==30,allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>%
    ggplot(aes(x=deletions,y=homeology_rate)) + geom_point() + scale_x_log10() + facet_wrap(~allelic_status_brca1_brca2,nrow=1)
dir.create(paste0(outDir,"/plots_homeology_rate_vs_deletions"))
ggsave(paste0(outDir,"/plots_homeology_rate_vs_deletions/homeology_rate_vs_deletions.pdf"),p,width=10,height=5)


## --- check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
make_distance_to_breakpoint_plot(candidates,outDir=outDir,similarities = c(80), alignment_lens = c(30))


## --- analyze repeats

# for candidates
candidates <- fread(paste0(outDir,"/ssa_events_candidates.tsv"))
candidates <- annotate_homeology_sequence_to_repeats(candidates)
candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.add_repeats.tsv"),sep="\t")

summarize_repeats(candidates,outDir=paste0(outDir,"/results_repeats"))
summarize_repeats_for_deletions(candidates %>% filter(min_similarity==80,alignment_length_cutoff==30),outDir=paste0(outDir,"/results_repeats/deletions_with_homeology"))

# for all deletions
deletions_add_repeats <- annotate_breakpoint_flanking_sequence_to_repeats(deletions)
deletions_add_repeats %>% fwrite(paste0(outDir,"/serena_BRCA-EU_SV_deletions.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))


## --- compare deletion length for deletions with/without homeology
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_30bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=30,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))


## -- check POLQ expression in HRDetect high and low groups
check_POLQ_expression_in_HRDetect_group(outDir)
check_POLQ_expression_in_HRDetect_group_excluding_BRCA1_BRCA2(outDir)

## -- check RAD52 expression in HRDetect high and low groups
check_RAD52_expression_in_HRDetect_group(outDir)
check_RAD52_expression_in_HRDetect_group_excluding_BRCA1_BRCA2(outDir)


##
## run for different deletion length
##

make_plots <- function(deletions, outDir, similarities = c(80), alignment_lens = c(25,30)){

    examine_deletion_length(
        deletions %>% 
            left_join(
                get_full_tumor_list() %>% 
                dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
                dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
        outDir = paste0(outDir,"/deletion_length"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))
        
    examine_deletion_length(
        deletions %>% 
            left_join(
                get_full_tumor_list() %>% 
                dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
                dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
        outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA2','control')))

    blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>%
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length)
        
    # get SSA homeology candidates by similarity + homeology length
    out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = similarities, alignment_lens = alignment_lens)

    candidates <- out[[1]]
    sample     <- out[[2]]

    candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
    sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

}


##
## for deletions 1-10kb 
##

outDir=paste0(outdir,"/output_deletions_1KBto10KB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1000, del_length<=10000)

make_plots(deletions,outDir)


##
## for deletions >10kb 
##

outDir=paste0(outdir,"/output_deletions_above_10KB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>10000)

make_plots(deletions,outDir)


##
## for deletions 10-100kb 
##

outDir=paste0(outdir,"/output_deletions_10KBto100KB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e4, del_length<=1e5)
nrow(deletions)

make_plots(deletions,outDir)


##
## for deletions >100kb 
##

outDir=paste0(outdir,"/output_deletions_above_100KB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e5)
nrow(deletions)

make_plots(deletions,outDir)


##
## for deletions 100kb-1Mb 
##

outDir=paste0(outdir,"/output_deletions_100KBto1MB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e5, del_length<=1e6)
nrow(deletions)

make_plots(deletions,outDir)


##
## for deletions 1Mb-10Mb
##

outDir=paste0(outdir,"/output_deletions_1MBto10MB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e6, del_length<=1e7)
nrow(deletions)

make_plots(deletions,outDir)

##
## for deletions >10Mb
##

outDir=paste0(outdir,"/output_deletions_above_10MB_mostSimilarAlignment_updatedControls"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
                dplyr::filter(del_length>1e7)
nrow(deletions)

make_plots(deletions,outDir)


##
## make homeology rate plot for deletion length 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb, keep y-axis in the same scale
##

outDir=paste0(outdir,"/output_deletions_above_1KB_mostSimilarAlignment_updatedControls/homeology_rate_in_different_deletion_size"); dir.create(outDir)

dat1 <- fread(paste0(outdir,"/output_deletions_1KBto10KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1-10kb')
dat2 <- fread(paste0(outdir,"/output_deletions_10KBto100KB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='10-100kb')
dat3 <- fread(paste0(outdir,"/output_deletions_100KBto1MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='100kb-1Mb')
dat4 <- fread(paste0(outdir,"/output_deletions_1MBto10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1Mb-10Mb')
dat5 <- fread(paste0(outdir,"/output_deletions_above_10MB_mostSimilarAlignment_updatedControls/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='>10Mb')

dat <- rbind(dat1,dat2,dat3,dat4,dat5) %>% dplyr::mutate(deletion_size=factor(deletion_size,levels=c('1-10kb','10-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')))
fwrite(dat,paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

# group by deletion size and genotype, BRCA1/BRCA2/control, 25bp
make_homeology_rate_plot_in_deletion_sizes(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==25),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_in_brca1_brca2_tumors.25bp.pdf",
    label.y=c(0.35,0.38,0.41), ylim=c(0,0.5))
# group by deletion size and genotype, BRCA2/control, 25bp
make_homeology_rate_plot_in_deletion_sizes(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==25),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_in_brca2_tumors.25bp.pdf",
    label.y=c(0.35,0.38), ylim=c(0,0.5))

# group by deletion size and genotype, BRCA1/BRCA2/control, 25bp
make_homeology_rate_plot_in_deletion_sizes(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_in_brca1_brca2_tumors.30bp.pdf",
    label.y=c(0.22,0.25,0.28), ylim=c(0,0.35))
# group by deletion size and genotype, BRCA2/control, 25bp
make_homeology_rate_plot_in_deletion_sizes(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_in_brca2_tumors.30bp.pdf",
    label.y=c(0.22,0.25), ylim=c(0,0.35))


# group by deletion size, BRCA1/BRCA2/control 25bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==25),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca1_brca2_tumors.25bp.pdf", ylim=c(0,0.5))
# group by deletion size, BRCA2/control 25bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==25),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2_tumors.25bp.pdf", ylim=c(0,0.5))
# group by deletion size, BRCA2 only 25bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==25),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2-only_tumors.25bp.pdf", ylim=c(0,0.5))

# group by deletion size, BRCA1/BRCA2/control 30bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca1_brca2_tumors.30bp.pdf", ylim=c(0,0.3))
# group by deletion size, BRCA2/control 30bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2','control')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2_tumors.30bp.pdf", ylim=c(0,0.3))
# group by deletion size, BRCA2 only 30bp
make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat %>% dplyr::filter(min_similarity==80, alignment_length_cutoff==30),outDir,group="allelic_status_brca1_brca2",
    comparisons = list(c('BRCA2')),
    outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2-only_tumors.30bp.pdf", ylim=c(0,0.3))


## -----------------------------------------------------------------------------------------------








## -----------------------------------------------------------------------------------------------
# below are the analysis not for the paper but for the further investigation


#
# ER+ tumors
#

outDir=paste0(outdir,"/output_deletions_above_1KB_mostSimilarAlignment_erpos_updatedControls"); dir.create(outDir)

sample_info <- get_full_tumor_list() %>% dplyr::filter(tumor_subtype=='ER+')

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
             dplyr::filter(del_length>1000, SAMPLE.TUMOR %in% sample_info$sampleid)

examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))
examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA2','control')))

blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>%
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length) %>%
    dplyr::filter(SAMPLE.TUMOR %in% sample_info$sampleid)

# get SSA homeology candidates by similarity + homeology length
out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))

candidates <- out[[1]]
sample     <- out[[2]] %>% dplyr::filter(tumor_subtype=='ER+')

candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

# for similarity x homeology seqeunce length
make_sample_plots_similarity_vs_alignment_length(sample=sample,outDir=outDir)

# for 80% similarity
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80"),
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

# for 80% similarity, BRCA2/control
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'),min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80_BRCA2_control"),
    comparisons=list(c('BRCA2','control')))


# check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
make_distance_to_breakpoint_plot(candidates,outDir=outDir,similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))


# analyze repeats
candidates <- fread(paste0(outDir,"/ssa_events_candidates.tsv"))
candidates <- annotate_homeology_sequence_to_repeats(candidates)
candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.add_repeats.tsv"),sep="\t")

summarize_repeats(candidates,outDir=paste0(outDir,"/results_repeats"))

# compare deletion length for deletions with/without homeology
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_30bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=30,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))
    
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_25bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=25,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))


#
# Triple negative tumors
#

outDir=paste0(outdir,"/output_deletions_above_1KB_mostSimilarAlignment_trineg"); dir.create(outDir)

sample_info <- get_full_tumor_list() %>% dplyr::filter(tumor_subtype=='Triple-')

deletions <- fread(paste0(outdir,"/serena_BRCA-EU_SV_deletions.tsv")) %>%
             dplyr::filter(del_length>1000, SAMPLE.TUMOR %in% sample_info$sampleid)

examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

examine_deletion_length(
    deletions %>% 
        left_join(
            get_full_tumor_list() %>% 
            dplyr::mutate(SAMPLE.TUMOR=sample_name) %>% 
            dplyr::select(SAMPLE.TUMOR,tumor_subtype,allelic_status_brca1_brca2)),
    outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
    group="allelic_status_brca1_brca2",
    comparisons=list(c('BRCA2','control')))

blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>% 
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length) %>%
    dplyr::filter(SAMPLE.TUMOR %in% sample_info$sampleid)

# get SSA homeology candidates by similarity + homeology length
out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))

candidates <- out[[1]]
sample     <- out[[2]] %>% dplyr::filter(tumor_subtype=='Triple-')

candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

# for similarity x homeology seqeunce length
make_sample_plots_similarity_vs_alignment_length(sample=sample,outDir=outDir)

# for 80% similarity
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80"),
    group = "allelic_status_brca1_brca2",
    comparisons = list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

# for 80% similarity, BRCA2/control
make_sample_plots_alignment_length(
    sample = sample %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'),min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30)),
    outDir=paste0(outDir,"/plots_similarity80_BRCA2_control"),
    comparisons=list(c('BRCA2','control')))
    

# check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
make_distance_to_breakpoint_plot(candidates,outDir=outDir,similarities = c(0,80,100), alignment_lens = c(10,20,25,30,40))


# analyze repeats
candidates <- fread(paste0(outDir,"/ssa_events_candidates.tsv"))
candidates <- annotate_homeology_sequence_to_repeats(candidates)
candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.add_repeats.tsv"),sep="\t")

summarize_repeats(candidates,outDir=paste0(outDir,"/results_repeats"))


# compare deletion length for deletions with/without homeology
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_30bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=30,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))
    
compare_deletion_length_for_deletions_with_and_without_homeology(
    deletions,
    candidates = fread(paste0(outDir,"/ssa_events_candidates.tsv")),
    outDir=paste0(outDir,"/deletion_length_for_deletions_with_and_without_homeology/results_similarity_80pct_25bp"),
    genotypes=c('BRCA2','control'),
    min_similarity=80, alignment_length_cutoff=25,
    comparisons = list(c('With homeology BRCA2','Without homeology BRCA2'),
                           c('With homeology BRCA2','With homeology control'),
                           c('Without homeology BRCA2','Without homeology control'),
                           c('With homeology control','Without homeology control')))