# Predict structrual variants (SV) with homeology in Serena's 560 breast cancer
rm(list = ls())
library(data.table)
library(dplyr)

# set output directory
outdir="example/input"; dir.create(outdir)


## -----------------------------------------------------------------------------------------------
# prepare SV deletions from Serena data

prepare_data <- function(){

    # Read in the data
    d <- fread("data/serena_data/BRCA-EU_SV/structural_somatic_mutation.BRCA-EU.tsv")
    table(d$variant_type)

    # Filter for deletions
    d <- d[variant_type == "deletion"]

    # select 5 random samples
    set.seed(1234)
    samples <- sample(unique(d$submitted_sample_id), 5)
    d <- d[submitted_sample_id %in% samples]

    # Encoding fields for the pipeline calculation
    d$sampleid <- gsub("(a|b).*$", "", d$submitted_sample_id)
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

    # write out the example data
    d.out %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position) %>%
        fwrite("example/input/structural_deletions_example.tsv",sep="\t")

}

prepare_data()
