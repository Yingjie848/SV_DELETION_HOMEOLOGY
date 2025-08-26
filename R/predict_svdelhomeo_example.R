# Predict structrual variants (SV) with homeology in Serena's 560 breast cancer
rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(tidyr)

source("lib/lib_blast.R")
source("lib/lib_identify_homeology.R")
source("lib/lib_plots.R")
source("lib/lib_repeats.R")
source("lib/lib_example_data.R")


## -----------------------------------------------------------------------------------------------
# Step 1: set up I/O files and directories

# set deletions input file
deletions <- fread("example/input/structural_deletions_example.tsv") %>%
            dplyr::mutate(del_length = end_position - start_position + 1) %>%
            dplyr::filter(del_length>1000)

# set output directory
outdir="example/output"; dir.create(outdir)

# set blast directory
blast_dir=paste0(outdir,"/blast_output"); dir.create(blast_dir)


## -----------------------------------------------------------------------------------------------
# Step 2: Run blastn, align 100bp leftside and 100bp rightside of the first breakpoint to 100bp 
# leftside and 100bp rightside of the second breakpoint

# run pairwise blast for two sequences of each deletion
# if everything runs well, delete tmp_dir afterwards, which contains intermediate blast output files
run_pairwise_blast_lr100(
    deletion_tbl=deletions,
    bltbl_file=paste0(blast_dir,"/blast_output.tsv"),
    tmp_dir=paste0(blast_dir,"/tmp"))

# clean up the folder
system(paste0("rm -rf ",blast_dir,"/tmp*"))

# Process blast output, extract alignments with the same direction
# Output: blast_output.same_direction.tsv
process_blast_output(blast_dir)


## -----------------------------------------------------------------------------------------------
# Step 3: Identify deletions with homeology candidates
#   Until this step, we have got the alignments, next we will identify deletions with homeology 
#   candidates by taking most similar and longest alignment, where the alignment should have >=80% 
#   similarity, and >=30bp alignment length

# load blast output data with the same direction, add deletion length
blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.same_direction.tsv")) %>%
    dplyr::left_join(
        deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length) %>% distinct, 
        by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
    dplyr::filter(!is.na(del_length)) %>%
    dplyr::select(-del_length)
    
# get SSA homeology candidates by similarity + homeology length
# For the final analysis, we only use 80% similarity and 30bp alignment length
similarities <- c(80)
alignment_lens <- c(30)
out <- get_candidates(
        full_deletions = deletions,
        blast_output_same_strand = blast_output_same_strand,
        outDir = outdir,
        take_longest = FALSE, 
        take_most_similar = TRUE, 
        similarities = similarities,
        alignment_lens = alignment_lens)

candidates <- out[[1]]
sample     <- out[[2]]

# save deletions with homeology
candidates %>% fwrite(paste0(outdir,"/ssa_events_candidates.tsv"),sep="\t")

# save summary of samples with homeology 
# (average values of deletion length, alignment length, identify/similarity, and homeology rate result, etc)
sample     %>% fwrite(paste0(outdir,"/ssa_events_samples.tsv"),sep="\t")


