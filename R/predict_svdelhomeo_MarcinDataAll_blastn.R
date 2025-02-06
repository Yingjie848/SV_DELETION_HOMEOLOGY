rm(list = ls())
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggforce)
library(ggpubr)

source("lib/lib_blast.R")
source("lib/lib_plots.R")
source("lib/lib_identify_homeology.R")
source("lib/lib_repeats.R")
source("lib/lib_marcin_data.R")


outdir="output/marcin_allTumor_lr100_blast/"; dir.create(outdir)

blast_dir=paste0(outdir,"/blast_output"); dir.create(blast_dir)


## -----------------------------------------------------------------------------------------------
# prepare deletions from Marcin's data

prepare_data <- function(){

    # Read junctions and hrd table data
    junctions.dt <- readRDS("input/setton_hadi_choo_2023/revision.junctions.dt.rds")
    hrd_tbl      <- readRDS("input/setton_hadi_choo_2023/hrd-supp-table_ZC.rds")

    # extract data
    hrd_tbl %>% as.data.frame %>% fwrite("input/setton_hadi_choo_2023/hrd-supp-table_ZC.tsv",sep="\t")
    hrd_tbl %>% as.data.frame %>% dplyr::filter(dataset=='PCAWG') %>% fwrite("input/setton_hadi_choo_2023/hrd-supp-table_ZC.PCAWG.tsv",sep="\t")
    hrd_tbl %>% as.data.frame %>% dplyr::filter(dataset=='TCGA') %>% fwrite("input/setton_hadi_choo_2023/hrd-supp-table_ZC.TCGA.tsv",sep="\t")
    hrd_tbl %>% as.data.frame %>% dplyr::filter(dataset=='PCAWG',fmut_bi %in% c('BRCA1','BRCA2','WT')) %>% fwrite("input/setton_hadi_choo_2023/hrd-supp-table_ZC.PCAWG_BRCA1_BRCA2_WT.tsv",sep="\t")

    hrd_tbl %>% filter(pair %in% junctions.dt$sample) %>% 
        group_by(dataset,fmut_bi) %>% 
        summarise(n=length(pair)) %>% as.data.frame %>% 
        fwrite(paste0(outdir,"/marcin_data_summary.tsv"),sep="\t")

    # check BRCA1/BRCA2/WT samples with deletions
    j <- junctions.dt %>% left_join(hrd_tbl %>% dplyr::mutate(sample=pair) %>% 
                        dplyr::select(sample,fmut_bi,dataset,HRDetect,biallelic_pathogenic_tier2,biallelic_likpathogenic_tier2),by='sample')
    j %>% filter(dataset %in% c('TCGA','PCAWG'),fmut_bi %in% c("BRCA1", "BRCA2", "WT"),!is.na(del)) %>% pull(sample) %>% unique %>% length
    j %>% filter(dataset %in% c('TCGA','PCAWG'),fmut_bi %in% c("BRCA1", "BRCA2", "WT"),!is.na(del),  span>1e3) %>% pull(sample) %>% unique %>% length

    # get deletions
    dels.dt = junctions.dt[(!is.na(del)) &
                        (span > 1e3) & 
                        (hrd_tbl[sample, fmut_bi %in% c("BRCA1", "BRCA2", "WT")])]
    dels.dt[, SAMPLE.TUMOR   := sample]
    dels.dt[, CHROM          := seqnames.1]
    dels.dt[, start_position := start.1]
    dels.dt[, end_position   := start.2]
    dels.dt[, del_length     := span]
    # add BRCA1/BRCA2/WT mutations
    dels.dt[, genotype := hrd_tbl[sample, fmut_bi]]
    dels.dt <- dels.dt %>% left_join(hrd_tbl %>% dplyr::mutate(sample=pair) %>% 
                        dplyr::select(sample,fmut_bi,dataset,HRDetect,biallelic_pathogenic_tier2,biallelic_likpathogenic_tier2),by='sample') %>%
                        dplyr::mutate('allelic_status_brca1_brca2'=gsub("WT","control",fmut_bi)) %>%
                        dplyr::rename('HRDetect_score'='HRDetect') %>%
                        dplyr::mutate(
                                HRDetect_group = ifelse(HRDetect_score < 0.2, "HRDetect-low",
                                    ifelse(HRDetect_score >= 0.7, "HRDetect-high", "HRDectect-intermediate")))

    summary(dels.dt$span)
    summary(as.numeric(dels.dt$hlen))
    table(dels.dt$allelic_status_brca1_brca2)


    # get deletions with homeology by samples predicted by Marcin's group, code copied from https://github.com/mskilab/setton_hadi_choo_2023/blob/main/notebooks/figures.ipynb
    plot_frac <- function(hlen_cutoff=30){
            hdels.dt = dels.dt[, .(ndels = .SD[, .N], nhdels = .SD[ hlen > hlen_cutoff, .N]), keyby = sample]
            hdels.dt[, genotype := hrd_tbl[sample, fmut_bi]]
            hdels.dt[, dataset := hrd_tbl[sample, dataset]]
            hdels.dt[, frac := nhdels / ndels]
            nrow(hdels.dt)
            table(hdels.dt$genotype)

            hdels.dt %>% group_by(genotype) %>% summarise(nhdels=sum(nhdels),mean_homeology_rate=mean(frac,na.rm=TRUE),median_homeology_rate=median(frac,na.rm=TRUE))

            # set levels for plotting
            hdels.dt[, genotype_level := ordered(genotype, levels = c("BRCA1", "BRCA2", "WT"))]
            hdels.dt[, dataset_level := ifelse(dataset == "STARR", "STARR", "OTHER")]

            # plot deletions with homeology by samples
            #color.map = c("STARR" = alpha("black", 0.6), "OTHER" = alpha("#a1a0a0", 0.6))

            pt = ggplot(hdels.dt, aes(x = genotype_level, y = frac, group = genotype_level)) + 
                    geom_violin(position = "identity", bw = 0.01) +
                    ggpubr::stat_compare_means(comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','WT'),c('BRCA2','WT')), 
                                    method="wilcox.test",method.args=list('exact'=FALSE)) +
                    ggforce::geom_sina(aes(x = genotype_level, y = frac, group = genotype), 
                                    position = "identity",
                                    size = 0.6) +
                    #scale_color_manual(values = color.map) +
                    ggpubr::theme_pubr() +
                    #ylim(0,0.3) +
                    labs(x = "",y = "fraction of deletions",title=hlen_cutoff) + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust=1))
    }

    plots <- lapply(c(0,10,20,30,40),plot_frac)
    p <- ggarrange(plotlist=plots,nrow=1,common.legend=TRUE)
    pdf(paste0(outdir,"/marcin_frac_hdels.pdf"),width=10,height=5);print(p);dev.off()


    # Write out the data
    d.out <- dels.dt %>% dplyr::filter(CHROM!='hs37d5') # hs37d5 comes from 1000 Genome Project
    d.out %>% fwrite(paste0(outdir,"/dels.tsv"), sep = "\t")


    # make deletion length distribution
    p <- ggplot(dels.dt,aes(x=span,y=after_stat(count/sum(count)))) + 
            geom_histogram(bins=50,fill="gray",color="blue") + 
            scale_x_log10(breaks=c(1,10,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)) + 
            xlab("Deletion length") + ylab("Frequency")
    ggsave(paste0(outdir,"/dels_histogram.pdf"),width=7,height=7)

    # make homeology length distribution
    p <- ggplot(dels.dt %>% dplyr::filter(hlen>0),aes(x=hlen,y=after_stat(count/sum(count)))) + 
            geom_histogram(bins=50,fill="gray",color="blue") + 
            #scale_x_log10(breaks=c(1,10,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)) + 
            xlab("Homeology length") + ylab("Frequency")
    ggsave(paste0(outdir,"/hdels_histogram.pdf"),width=7,height=7)

    # count number of deletions per sample
    deletions_per_sample <- dels.dt %>% group_by(sample) %>% dplyr::summarise(deletions=length(sample))
    deletions_per_sample %>% fwrite(paste0(outdir,"/dels.dels_per_sample.tsv"),sep="\t")

    d.out

}

# get deletions from tumors with BRCA1/2/WT mutations
d.out <- prepare_data() 


## -----------------------------------------------------------------------------------------------
# Run blastn, aligning 100bp leftside and 100bp rightside of the first breakpoint to  
#   100bp leftside and 100bp rightside of the second breakpoint

# run pairwise blast for two sequences of each deletion
# if everything runs well, delete tmp_dir afterwards, which contains intermediate blast output files
run_pairwise_blast_lr100(
    deletion_tbl=d.out,
    bltbl_file=paste0(blast_dir,"/blast_output.tsv"),
    tmp_dir=paste0(blast_dir,"/tmp"))

# Process blast output, extract alignments with the same direction
# Output: blast_output.same_direction.tsv
process_blast_output_marcin_data(blast_dir)


## -----------------------------------------------------------------------------------------------
# To this step, we have made the alignment, next we identify deletions with homeology candidates 
#   by taking most similar and longest alignment, where the alignments should have >=80% similarity, 
#   and >=30bp alignment length

make_ssa_candidate_table <- function(deletions, outDir, similarities = c(80), alignment_lens = c(30)){

    # examine deletion length for BRCA1/2/control tumors
    examine_deletion_length(
        deletions,
        outDir = paste0(outDir,"/deletion_length"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA1','BRCA2'),c('BRCA1','control'),c('BRCA2','control')))

    # examine deletion length for BRCA2/control tumors
    examine_deletion_length(
        deletions,
        outDir = paste0(outDir,"/deletion_length_BRCA2_control"),
        group="allelic_status_brca1_brca2",
        comparisons=list(c('BRCA2','control')))

    # load blast output data with the same direction, add deletion length, and filter out deletions without deletion length
    blast_output_same_strand <- fread(paste0(blast_dir,"/blast_output.filtered_same_direction.tsv")) %>%
                                dplyr::left_join(deletions %>% dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position,del_length), 
                                                by=c('SAMPLE.TUMOR','CHROM','start_position','end_position')) %>%
                                dplyr::filter(!is.na(del_length)) %>%
                                dplyr::select(-del_length)
        
    # get SSA homeology candidates by similarity + homeology length
    # For the final analysis, we only use 80% similarity and 30bp alignment length
    out <- get_candidates(deletions,blast_output_same_strand,outDir,take_longest = FALSE, take_most_similar = TRUE, similarities = similarities, alignment_lens = alignment_lens)

    candidates <- out[[1]] # deletions with homeology candidates
    sample     <- out[[2]] # sample level results, including homeology rate (proportion of deletions with homeology)

    candidates %>% fwrite(paste0(outDir,"/ssa_events_candidates.tsv"),sep="\t")
    sample     %>% fwrite(paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

    out

}


##
## for deletions >=1kb
##

outDir=paste0(outdir,"/output_deletions_above_1KB"); dir.create(outDir)

## --- extract deletions
deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
     dplyr::filter(del_length>=1e3)
nrow(deletions) # 47336

## --- count deletions in BRCA2 and control
dels <- deletions %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 45319

# check number of deletions in a range of sizes: 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
dels %>% filter(del_length>1e3,del_length<=1e4) %>% nrow # 11248
dels %>% filter(del_length>1e4,del_length<=1e5) %>% nrow # 15356
dels %>% filter(del_length>1e5,del_length<=1e6) %>% nrow # 13195
dels %>% filter(del_length>1e6,del_length<=1e7) %>% nrow # 3397
dels %>% filter(del_length>1e7) %>% nrow # 2123

## --- make SSA candidate table
out <- make_ssa_candidate_table(deletions,outDir)
candidates <- out[[1]]
sample     <- out[[2]]


## --- further analysis

## --- Grid search similarity x homeology seqeunce length
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


## --- check controls with long tails, it could be caused by low number of deletions
p <- sample %>% dplyr::filter(min_similarity==80,alignment_length_cutoff==30,allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>%
    ggplot(aes(x=deletions,y=homeology_rate)) + geom_point() + scale_x_log10() + facet_wrap(~allelic_status_brca1_brca2,nrow=1)
dir.create(paste0(outDir,"/plots_homeology_rate_vs_deletions"))
ggsave(paste0(outDir,"/plots_homeology_rate_vs_deletions/homeology_rate_vs_deletions.pdf"),p,width=10,height=5)


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
deletions_add_repeats %>% fwrite(paste0(outDir,"/dels.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))

## --- make boxplot for different tumor type
sample <- fread(paste0(outDir,"/ssa_events_samples.tsv"))
make_plots_for_tumor_types_allTumors(sample,hrd_tbl,outDir=paste0(outDir,"/plots_tumor_type_homeology_rate"))
make_plots_for_tumor_types_allBOPP(sample,hrd_tbl,outDir=paste0(outDir,"/plots_tumor_type_homeology_rate_allBOPP"))



##
## for deletions 1KB-100KB
##

outDir=paste0(outdir,"/output_deletions_1KB-100KB_final"); dir.create(outDir)

## --- extract deletions
deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
     dplyr::filter(del_length>1e3, del_length<1e5)
nrow(deletions) # 27946

## --- count deletions in BRCA2 and control
dels <- deletions %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 26604

## --- make SSA candidate table
out <- make_ssa_candidate_table(deletions,outDir)
candidates <- out[[1]]
sample     <- out[[2]]


## --- further analysis

## --- Grid search similarity x homeology seqeunce length
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
deletions_add_repeats %>% fwrite(paste0(outDir,"/dels.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))

## --- make boxplot for different tumor type
sample <- fread(paste0(outDir,"/ssa_events_samples.tsv"))
make_plots_for_tumor_types_allTumors(sample,hrd_tbl,outDir=paste0(outDir,"/plots_tumor_type_homeology_rate"))
make_plots_for_tumor_types_allBOPP(sample,hrd_tbl,outDir=paste0(outDir,"/plots_tumor_type_homeology_rate_allBOPP"))



##
## run for different deletion length
## 


##
## for deletions 1-10kb
##

outDir=paste0(outdir,"/output_deletions_1KBto10KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
     dplyr::filter(del_length>=1e3, del_length<=1e4)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions >10kb 
##

outDir=paste0(outdir,"/output_deletions_above_10KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>10000)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 10-100kb 
##

outDir=paste0(outdir,"/output_deletions_10KBto100KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>1e4, del_length<=1e5)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions >100kb 
##

outDir=paste0(outdir,"/output_deletions_above_100KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>1e5)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 100kb-1Mb 
##

outDir=paste0(outdir,"/output_deletions_100KBto1MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>1e5, del_length<=1e6)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 1Mb-10Mb
##

outDir=paste0(outdir,"/output_deletions_1MBto10MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>1e6, del_length<=1e7)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)

##
## for deletions >10Mb
##

outDir=paste0(outdir,"/output_deletions_above_10MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
                dplyr::filter(del_length>1e7)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)



##
## make homeology rate plot for deletion length 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb, keep y-axis in the same scale
##

make_homeology_rate_plot_for_different_deletion_sizes <- function(){

    outDir=paste0(outdir,"/output_deletions_above_1KB/homeology_rate_in_different_deletion_size"); dir.create(outDir)

    dat1 <- fread(paste0(outdir,"/output_deletions_1KBto10KB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1-10kb')
    dat2 <- fread(paste0(outdir,"/output_deletions_10KBto100KB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='10-100kb')
    dat3 <- fread(paste0(outdir,"/output_deletions_100KBto1MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='100kb-1Mb')
    dat4 <- fread(paste0(outdir,"/output_deletions_1MBto10MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1Mb-10Mb')
    dat5 <- fread(paste0(outdir,"/output_deletions_above_10MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='>10Mb')

    dat <- rbind(dat1,dat2,dat3,dat4,dat5) %>% dplyr::mutate(deletion_size=factor(deletion_size,levels=c('1-10kb','10-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')))
    fwrite(dat,paste0(outDir,"/ssa_events_samples.tsv"),sep="\t")

    # group by deletion size and genotype
    make_homeology_rate_plot_in_deletion_sizes(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_in_brca1_brca2_tumors.pdf",
        label.y=c(0.22,0.25,0.28), ylim=c(0,0.35))


    make_homeology_rate_plot_in_deletion_sizes(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_in_brca2_tumors.pdf",
        label.y=c(0.22,0.25), ylim=c(0,0.3))

    # group by deletion size
    make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca1_brca2_tumors.pdf", ylim=c(0,0.3))

    make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2_tumors.pdf", ylim=c(0,0.3))

    make_homeology_rate_plot_in_deletion_sizes_on_x_axis(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2')),
        outfile="homeology_rate_plot_in_deletion_sizes_on_x_axis_in_brca2-only_tumors.pdf", ylim=c(0,0.3))

}

make_homeology_rate_plot_for_different_deletion_sizes()


######################################################################################################
# Identify deletions with homeology only for BOPP tumors: BRCA, OV, PACA, PRAD
######################################################################################################

##
## for deletions >=1kb
##

outDir=paste0(outdir,"/BOPP_deletions_above_1KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>=1e3)
nrow(deletions) # 13771

## --- count deletions in BRCA2 and control
dels <- deletions %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 12178

# group1
dels %>% filter(del_length>1e3, del_length<=1e4) %>% nrow # 3481
dels %>% filter(del_length>1e4) %>% nrow # 8697

# group2 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb
dels %>% filter(del_length>1e3,del_length<=1e4) %>% nrow # 3481
dels %>% filter(del_length>1e4,del_length<=1e5) %>% nrow # 3586
dels %>% filter(del_length>1e5,del_length<=1e6) %>% nrow # 2810
dels %>% filter(del_length>1e6,del_length<=1e7) %>% nrow # 1413
dels %>% filter(del_length>1e7) %>% nrow # 888

# group3
dels %>% filter(del_length>1e5) %>% nrow # 5111

## --- make SSA candidate table
out <- make_ssa_candidate_table(deletions,outDir)
candidates <- out[[1]]
sample     <- out[[2]]

## --- for similarity x homeology seqeunce length
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


## --- check controls with long tails, it could be caused by low number of deletions
p <- sample %>% dplyr::filter(min_similarity==80,alignment_length_cutoff==30,allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>%
    ggplot(aes(x=deletions,y=homeology_rate)) + geom_point() + scale_x_log10() + facet_wrap(~allelic_status_brca1_brca2,nrow=1)
dir.create(paste0(outDir,"/plots_homeology_rate_vs_deletions"))
ggsave(paste0(outDir,"/plots_homeology_rate_vs_deletions/homeology_rate_vs_deletions.pdf"),p,width=10,height=5)


# check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
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
deletions_add_repeats %>% fwrite(paste0(outDir,"/dels.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))



##
## for deletions 1KB-100KB
##

outDir=paste0(outdir,"/BOPP_deletions_1KB-100KB_final"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>1e3, del_length<1e5)
nrow(deletions) # 8147

## --- count deletions in BRCA2 and control
dels <- deletions %>% dplyr::filter(allelic_status_brca1_brca2 %in% c('BRCA2','control'))
nrow(dels) # 7067

## --- make SSA candidate table
out <- make_ssa_candidate_table(deletions,outDir)
candidates <- out[[1]]
sample     <- out[[2]]

## --- for similarity x homeology seqeunce length
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


## --- check controls with long tails, it could be caused by low number of deletions
p <- sample %>% dplyr::filter(min_similarity==80,alignment_length_cutoff==30,allelic_status_brca1_brca2 %in% c('BRCA1','BRCA2','control')) %>%
    ggplot(aes(x=deletions,y=homeology_rate)) + geom_point() + scale_x_log10() + facet_wrap(~allelic_status_brca1_brca2,nrow=1)
dir.create(paste0(outDir,"/plots_homeology_rate_vs_deletions"))
ggsave(paste0(outDir,"/plots_homeology_rate_vs_deletions/homeology_rate_vs_deletions.pdf"),p,width=10,height=5)


# check homeology distance to breakpoints (distance means the nearest end of homeology sequence to breakpoint)
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
deletions_add_repeats %>% fwrite(paste0(outDir,"/dels.add_repeats.tsv"),sep="\t")
summarize_repeats_for_deletions(deletions_add_repeats,outDir=paste0(outDir,"/results_repeats/all_deletions"))

# for repeats in the whole genome
summarize_repeats_in_genome(outDir=paste0(outDir,"/results_repeats/genome_repeats"))

## --- make boxplot for different tumor type
sample <- fread(paste0(outDir,"/ssa_events_samples.tsv"))
make_plots_for_tumor_types_allBOPP(sample,hrd_tbl,outDir=paste0(outDir,"/plots_tumor_type_homeology_rate_allBOPP"))


## --- run for different deletion length

##
## for deletions 1-10kb
##

outDir=paste0(outdir,"/BOPP_deletions_1KBto10KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>1e3, del_length<1e4)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 10-100kb 
##

outDir=paste0(outdir,"/BOPP_deletions_10KBto100KB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>=1e4, del_length<1e5)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 100kb-1Mb 
##

outDir=paste0(outdir,"/BOPP_deletions_100KBto1MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>=1e5, del_length<1e6)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)


##
## for deletions 1Mb-10Mb
##

outDir=paste0(outdir,"/BOPP_deletions_1MBto10MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>=1e6, del_length<1e7)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)

##
## for deletions >10Mb
##

outDir=paste0(outdir,"/BOPP_deletions_above_10MB"); dir.create(outDir)

deletions <- fread(paste0(outdir,"/dels.tsv")) %>%
    left_join(hrd_tbl %>% dplyr::select(pair,tumor_type_final), by=c('SAMPLE.TUMOR'='pair')) %>%
    dplyr::filter(tumor_type_final %in% c('BRCA','OV','PACA','PRAD'), del_length>=1e7)
nrow(deletions)

make_ssa_candidate_table(deletions,outDir)



##
## make homeology rate plot for deletion length 1-10kb, 10-100kb, 100kb-1Mb, 1Mb-10Mb, >10Mb, keep y-axis in the same scale
##

make_homeology_rate_plot_for_different_deletion_sizes <- function(){

    outDir=paste0(outdir,"/BOPP_deletions_above_1KB/homeology_rate_in_different_deletion_size"); dir.create(outDir, recursive = TRUE)

    dat1 <- fread(paste0(outdir,"/BOPP_deletions_1KBto10KB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1-10kb')
    dat2 <- fread(paste0(outdir,"/BOPP_deletions_10KBto100KB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='10-100kb')
    dat3 <- fread(paste0(outdir,"/BOPP_deletions_100KBto1MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='100kb-1Mb')
    dat4 <- fread(paste0(outdir,"/BOPP_deletions_1MBto10MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='1Mb-10Mb')
    dat5 <- fread(paste0(outdir,"/BOPP_deletions_above_10MB/ssa_events_samples.tsv")) %>% dplyr::mutate(deletion_size='>10Mb')

    dat <- rbind(dat1,dat2,dat3,dat4,dat5) %>% dplyr::mutate(deletion_size=factor(deletion_size,levels=c('1-10kb','10-100kb','100kb-1Mb','1Mb-10Mb','>10Mb')))
    fwrite(dat,paste0(outDir,"/sample_homeology_rate.tsv"),sep="\t")

    # group by deletion size and genotype
    make_homeology_rate_plot_in_deletion_sizes(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control'),c('BRCA2','BRCA1'),c('BRCA1','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_in_brca1_brca2_tumors.pdf",
        label.y=c(0.22,0.25,0.28), ylim=c(0,0.35))


    make_homeology_rate_plot_in_deletion_sizes(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control')),
        outfile="homeology_rate_plot_in_deletion_sizes_in_brca2_tumors.pdf",
        label.y=c(0.22,0.25), ylim=c(0,0.3))


    make_homeology_rate_violinplot_in_deletion_sizes(dat,outDir,group="allelic_status_brca1_brca2",
        comparisons = list(c('BRCA2','control')),
        outfile="homeology_rate_violinplot_in_deletion_sizes_in_brca2_tumors.pdf",
        label.y=c(0.22,0.25), ylim=c(0,0.3))

}

make_homeology_rate_plot_for_different_deletion_sizes()