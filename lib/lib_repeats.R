## Input predicted deletions with homeology
## Output same table adding repeat information
## Input columns:
#   $ CHROM                     : chr  "1" "1" "1" "1" ...
#   $ start_position            : int  31536066 66759817 153165956 159254478 161070799 175351423 175507765 177070588 8790849 24561215 ...
#   $ end_position              : int  31844037 66770403 171834134 180618230 161120751 175541419 178981254 178446678 79365539 25051389 ...
## Output columns:
#   $ query_start               : num  3.15e+07 6.68e+07 1.53e+08 1.59e+08 1.61e+08 ...
#   $ query_end                 : num  3.15e+07 6.68e+07 1.53e+08 1.59e+08 1.61e+08 ...
#   $ subject_start             : num  3.18e+07 6.68e+07 1.72e+08 1.81e+08 1.61e+08 ...
#   $ subject_end               : num  3.18e+07 6.68e+07 1.72e+08 1.81e+08 1.61e+08 ...
#   $ query_rmsk_start          : int  NA 66759606 NA 159254331 NA 175350872 NA 177070263 NA 24561188 ...
#   $ query_rmsk_end            : int  NA 66759909 NA 159254584 NA 175351410 NA 177070833 NA 24561372 ...
#   $ query_repName             : chr  NA "AluSx1" NA "L1ME4a" ...
#   $ query_repClass            : chr  NA "SINE" NA "LINE" ...
#   $ query_repFamily           : chr  NA "Alu" NA "L1" ...
#   $ query_rmsk_overlap_width  : num  0 4 0 4 0 4 0 4 0 4 ...
#   $ subject_rmsk_start        : int  NA NA NA NA 161120434 NA 178981066 NA NA NA ...
#   $ subject_rmsk_end          : int  NA NA NA NA 161120741 NA 178981859 NA NA NA ...
#   $ subject_repName           : chr  NA NA NA NA ...
#   $ subject_repClass          : chr  NA NA NA NA ...
#   $ subject_repFamily         : chr  NA NA NA NA ...
#   $ subject_rmsk_overlap_width: num  0 0 0 0 4 0 4 0 0 0 ...
#   $ repeat_group              : chr  "No repeats" "One breakpoint" "No repeats" "One breakpoint" ...
## Notes: if you have query_start, query_end (1st breakpoint region), subject_start, subject_end (2nd breakpoint region), you don't need to calculate them in the function
annotate_homeology_sequence_to_repeats <- function(candidates){

    library(GenomicRanges)

    # intersect query region with subject region
    # for query, intersected with multiple subjects, take longest one
    get_overlap <- function(query_gr,subject_gr){
        
        fo   <- findOverlaps(query_gr, subject_gr)
        overlap_width <- width(pintersect(query_gr[queryHits(fo),],subject_gr[subjectHits(fo),]))
        ovl  <- as.data.frame(
                    cbind(data.frame(CHROM=as.vector(seqnames(query_gr))[queryHits(fo)],
                                        start_position=start(query_gr)[queryHits(fo)],
                                        end_position=end(query_gr)[queryHits(fo)]),
                                        sbj_start=start(subject_gr)[subjectHits(fo)],
                                        sbj_end=end(subject_gr)[subjectHits(fo)],
                          mcols(subject_gr)[subjectHits(fo),]
                    )) %>% 
                    dplyr::mutate(overlap_width=overlap_width) %>%
                    group_by(CHROM,start_position,end_position) %>% 
                    arrange(desc(overlap_width)) %>% 
                    slice_head(n=1)
        
    }

    # load repeatmasker and save to GRange object
    rmsk <- fread("data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz") %>% # /lila/data/riazlab/projects/zhuy1/my_projects/Manisha_SSA/ssa_homeology_20231119/data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz
        dplyr::select(genoName,genoStart,genoEnd,strand,repName,repClass,repFamily)
    rmsk_gr <- makeGRangesFromDataFrame(rmsk %>% dplyr::mutate(genoName=gsub("chr","",genoName)), 
                seqnames.field = 'genoName', start.field='genoStart', end.field='genoEnd',keep.extra.columns=T)
    rm(rmsk)

    # calculate query region start/end (1st breakpoint), subject region start/end (2nd breakpoint)
    # !!! you don't need to do this if you already have query_start, query_end, subject_start, subject_end
    candidates <- candidates %>% 
        dplyr::mutate(
            query_start=start_position-99+q_start-1,
            query_end  =start_position-99+q_end-1,
            subject_start=end_position-100+s_start-1,
            subject_end  =end_position-100+s_end-1
            )

    query_df   <- candidates %>% dplyr::select(CHROM,query_start,query_end) 
    subject_df <- candidates %>% dplyr::select(CHROM,subject_start,subject_end)

    # save query region (1st breakpoint sequence) to GRange object
    query_gr  <- makeGRangesFromDataFrame(query_df %>% distinct, 
                seqnames.field = 'CHROM', start.field='query_start', end.field='query_end',keep.extra.columns=T)
    
    # annotate query region to repeat region
    query_overlap <- get_overlap(query_gr,rmsk_gr)

    # add repeat information to query
    query_df <- query_df %>% 
        left_join(query_overlap %>% dplyr::rename(query_start=start_position,query_end=end_position),by=c('CHROM','query_start','query_end')) %>% 
        dplyr::mutate(overlap_width = ifelse(is.na(overlap_width),0,overlap_width)) %>%
        dplyr::rename(query_rmsk_start=sbj_start,query_rmsk_end=sbj_end,query_repName=repName,query_repClass=repClass,
                      query_repFamily=repFamily,query_rmsk_overlap_width=overlap_width)
    
    # do same thing for the second breakpoint sequence
    # map subject sequence to repeatmasker
    subject_gr  <- makeGRangesFromDataFrame(subject_df %>% distinct, 
                seqnames.field = 'CHROM', start.field='subject_start', end.field='subject_end',keep.extra.columns=T)
    
    subject_overlap <- get_overlap(subject_gr,rmsk_gr)

    subject_df <- subject_df %>% 
        left_join(subject_overlap %>% dplyr::rename(subject_start=start_position,subject_end=end_position),by=c('CHROM','subject_start','subject_end')) %>% 
        dplyr::mutate(overlap_width = ifelse(is.na(overlap_width),0,overlap_width)) %>%
        dplyr::rename(subject_rmsk_start=sbj_start,subject_rmsk_end=sbj_end,subject_repName=repName,subject_repClass=repClass,
                      subject_repFamily=repFamily,subject_rmsk_overlap_width=overlap_width)
    
    # add results to candidates table
    candidates <- cbind(candidates, query_df[,4:ncol(query_df)])
    candidates <- cbind(candidates, subject_df[,4:ncol(subject_df)])

    # group deletions by two breakpoints whether they overlapped with repeats
    candidates <- candidates %>% dplyr::mutate(
        repeat_group = ifelse(query_rmsk_overlap_width>0 & subject_rmsk_overlap_width>0, 'Both breakpoints',
                       ifelse(query_rmsk_overlap_width>0 | subject_rmsk_overlap_width>0, 'One breakpoint','No repeats')
    ))

    candidates

}

# annotate breakpoint flanking sequence to repeats  (100bp flanking sequence)
annotate_breakpoint_flanking_sequence_to_repeats <- function(deletions){

    library(GenomicRanges)

    # intersect query region with subject region
    # for query, intersected with multiple subjects, take longest one
    get_overlap <- function(query_gr,subject_gr){
        
        fo   <- findOverlaps(query_gr, subject_gr)
        overlap_width <- width(pintersect(query_gr[queryHits(fo),],subject_gr[subjectHits(fo),]))
        ovl  <- as.data.frame(
                    cbind(data.frame(CHROM=as.vector(seqnames(query_gr))[queryHits(fo)],
                                        start_position=start(query_gr)[queryHits(fo)],
                                        end_position=end(query_gr)[queryHits(fo)]),
                                        sbj_start=start(subject_gr)[subjectHits(fo)],
                                        sbj_end=end(subject_gr)[subjectHits(fo)],
                          mcols(subject_gr)[subjectHits(fo),]
                    )) %>% 
                    dplyr::mutate(overlap_width=overlap_width) %>%
                    group_by(CHROM,start_position,end_position) %>% 
                    arrange(desc(overlap_width)) %>% 
                    slice_head(n=1)
        
    }

    # load repeatmasker and save to GRange object
    rmsk <- fread("data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz") %>% # /lila/data/riazlab/projects/zhuy1/my_projects/Manisha_SSA/ssa_homeology_20231119/data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz
        dplyr::select(genoName,genoStart,genoEnd,strand,repName,repClass,repFamily)
    rmsk_gr <- makeGRangesFromDataFrame(rmsk %>% dplyr::mutate(genoName=gsub("chr","",genoName)), 
                seqnames.field = 'genoName', start.field='genoStart', end.field='genoEnd',keep.extra.columns=T)
    rm(rmsk)

    # calculate query region start/end (1st breakpoint), subject region start/end (2nd breakpoint)
    # !!! you don't need to do this if you already have query_start, query_end, subject_start, subject_end
    deletions <- deletions %>% 
        dplyr::mutate(
            query_start=start_position-99, # left 100bp flanking sequence, include first breakpoint position
            query_end  =start_position+100, # right 100bp flanking sequence, do not include first breakpoint position
            subject_start=end_position-100, # left 100bp flanking sequence, do not include second breakpoint position
            subject_end  =end_position+99 # right 100bp flanking sequence, include second breakpoint position
            )

    query_df   <- deletions %>% dplyr::select(CHROM,query_start,query_end) 
    subject_df <- deletions %>% dplyr::select(CHROM,subject_start,subject_end)

    # save query region (1st breakpoint sequence) to GRange object
    query_gr  <- makeGRangesFromDataFrame(query_df %>% distinct, 
                seqnames.field = 'CHROM', start.field='query_start', end.field='query_end',keep.extra.columns=T)
    
    # annotate query region to repeat region
    query_overlap <- get_overlap(query_gr,rmsk_gr)

    # add repeat information to query
    query_df <- query_df %>% 
        left_join(query_overlap %>% dplyr::rename(query_start=start_position,query_end=end_position),by=c('CHROM','query_start','query_end')) %>% 
        dplyr::mutate(overlap_width = ifelse(is.na(overlap_width),0,overlap_width)) %>%
        dplyr::rename(query_rmsk_start=sbj_start,query_rmsk_end=sbj_end,query_repName=repName,query_repClass=repClass,
                      query_repFamily=repFamily,query_rmsk_overlap_width=overlap_width)
    
    # do same thing for the second breakpoint sequence
    # map subject sequence to repeatmasker
    subject_gr  <- makeGRangesFromDataFrame(subject_df %>% distinct, 
                seqnames.field = 'CHROM', start.field='subject_start', end.field='subject_end',keep.extra.columns=T)
    
    subject_overlap <- get_overlap(subject_gr,rmsk_gr)

    subject_df <- subject_df %>% 
        left_join(subject_overlap %>% dplyr::rename(subject_start=start_position,subject_end=end_position),by=c('CHROM','subject_start','subject_end')) %>% 
        dplyr::mutate(overlap_width = ifelse(is.na(overlap_width),0,overlap_width)) %>%
        dplyr::rename(subject_rmsk_start=sbj_start,subject_rmsk_end=sbj_end,subject_repName=repName,subject_repClass=repClass,
                      subject_repFamily=repFamily,subject_rmsk_overlap_width=overlap_width)
    
    # add results to deletions table
    deletions <- cbind(deletions, query_df[,4:ncol(query_df)])
    deletions <- cbind(deletions, subject_df[,4:ncol(subject_df)])

    # group deletions by two breakpoints whether they overlapped with repeats
    deletions <- deletions %>% dplyr::mutate(
        repeat_group = ifelse(query_rmsk_overlap_width>0 & subject_rmsk_overlap_width>0, 'Both breakpoints',
                       ifelse(query_rmsk_overlap_width>0 | subject_rmsk_overlap_width>0, 'One breakpoint','No repeats')
    ))

    deletions

}


group_deletions_by_repeats <- function(candidates){
    # group deletions by two breakpoints whether they overlapped with repeats
    candidates <- candidates %>% dplyr::mutate(
        # consider repeats not LINE/SINE/LTR/DNA as others, here set query_rmsk_overlap_width and subject_rmsk_overlap_width to 0
        query_rmsk_overlap_width = ifelse(is.na(query_repClass) | !query_repClass %in% c('LINE','SINE','LTR','DNA'),0,query_rmsk_overlap_width),
        subject_rmsk_overlap_width = ifelse(is.na(subject_repClass) | !subject_repClass %in% c('LINE','SINE','LTR','DNA'),0,subject_rmsk_overlap_width),
        repeat_group = ifelse(query_rmsk_overlap_width>0 & subject_rmsk_overlap_width>0, 'Both breakpoints',
                       ifelse(query_rmsk_overlap_width>0 | subject_rmsk_overlap_width>0, 'One breakpoint','No repeats')
    ))
}

summarize_repeats <- function(candidates,outDir){

    library(tidyr)

    dir.create(outDir)

    # group deletions by two breakpoints whether they overlapped with repeats
    candidates <- group_deletions_by_repeats(candidates)

    repeats_summary <- candidates %>% 
        group_by(min_similarity,alignment_length_cutoff,allelic_status_brca1_brca2) %>% 
        summarise(Both=sum(repeat_group=="Both breakpoints"),
                One=sum(repeat_group=="One breakpoint"),
                No=sum(repeat_group=="No repeats"),
                total=Both+One+No,
                query_overlap_100pct     =sum(query_rmsk_overlap_width/(query_end-query_start+1)==1),
                query_overlap_above_90pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.9),
                query_overlap_above_80pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.8),
                query_overlap_above_50pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.5),
                subject_overlap_100pct     =sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)==1),
                subject_overlap_above_90pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.9),
                subject_overlap_above_80pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.8),
                subject_overlap_above_50pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.5)
                ) %>%
        dplyr::mutate(pct_Both=Both/total*100,pct_One=One/total*100,pct_No=No/total*100) %>%
        dplyr::mutate(
            pct_query_overlap_100pct     =query_overlap_100pct/(Both+One)*100,
            pct_query_overlap_above_90pct=query_overlap_above_90pct/(Both+One)*100,
            pct_query_overlap_above_80pct=query_overlap_above_80pct/(Both+One)*100,
            pct_query_overlap_above_50pct=query_overlap_above_50pct/(Both+One)*100,
            pct_subject_overlap_100pct     =subject_overlap_100pct/(Both+One)*100,
            pct_subject_overlap_above_90pct=subject_overlap_above_90pct/(Both+One)*100,
            pct_subject_overlap_above_80pct=subject_overlap_above_80pct/(Both+One)*100,
            pct_subject_overlap_above_50pct=subject_overlap_above_50pct/(Both+One)*100
        )
    repeats_summary  %>% fwrite(paste0(outDir,"/repeats_summary.tsv"),sep="\t")

    # make barplot for homeology length 10,20,25,30bp
    p <- repeats_summary %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(10,20,25,30), allelic_status_brca1_brca2 != 'other') %>%
        dplyr::mutate(group=paste0(min_similarity,'% ',alignment_length_cutoff,'bp')) %>%
        pivot_longer(cols=c('pct_Both','pct_One','pct_No'),names_to="repeats_location",values_to="pct_repeats") %>%
        dplyr::mutate(repeats_location=case_when(
            repeats_location=="pct_Both" ~ "Both breakpoints",
            repeats_location=="pct_One" ~ "One breakpoint",
            repeats_location=="pct_No" ~ "No repeats",
        )) %>%
        dplyr::mutate(repeats_location=factor(repeats_location,levels=c("Both breakpoints","One breakpoint","No repeats"))) %>%
        ggplot(aes(x=allelic_status_brca1_brca2,y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", position="fill") + 
            theme_bw() + 
            xlab('') + ylab('Proportion') + 
            ggtitle("Deletions around repeats (both/one/no breakpoints)") + 
            facet_wrap(~group,nrow=1)
    ggsave(paste0(outDir,"/repeats_summary.barplot.pdf"),width=10,height=5)


    # make barplot for homeology length 30bp

    repeats_summary_30bp <- repeats_summary %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(30), allelic_status_brca1_brca2 != 'other') %>%
        dplyr::mutate(group=paste0(min_similarity,'% ',alignment_length_cutoff,'bp')) %>%
        pivot_longer(cols=c('pct_Both','pct_One','pct_No'),names_to="repeats_location",values_to="pct_repeats") %>%
        dplyr::mutate(repeats_location=case_when(
            repeats_location=="pct_Both" ~ "Both breakpoints",
            repeats_location=="pct_One" ~ "One breakpoint",
            repeats_location=="pct_No" ~ "No repeats",
        )) %>%
        dplyr::mutate(repeats_location=factor(repeats_location,levels=c("Both breakpoints","One breakpoint","No repeats")))

    repeats_summary_30bp %>% fwrite(paste0(outDir,"/repeats_summary.homelen_30.tsv"),sep="\t")

    p <- repeats_summary_30bp %>%
        ggplot(aes(x=allelic_status_brca1_brca2,y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", position="fill") + 
            theme_bw() + 
            xlab('') + ylab('Proportion') + 
            ggtitle("Deletions around repeats (both/one/no breakpoints)")
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_barplot.pdf"),width=5,height=5)

    # make piechart using ggplot2
    dplot <- repeats_summary_30bp %>% group_by(min_similarity,alignment_length_cutoff,group,repeats_location) %>%
        summarise(pct_repeats=mean(pct_repeats),) %>% ungroup()
    dplot %>% fwrite(paste0(outDir,"/repeats_summary.homelen_30_piechart.tsv"),sep="\t")
    p <- dplot %>%
        ggplot(aes(x="",y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", width=1) + 
            coord_polar("y", start=0) + 
            geom_text(aes(label = paste0(round(pct_repeats,0),"%")), size = 6, color = "black",
                position = position_stack(vjust = 0.5)) +
            theme_void() + scale_fill_hue(name="Repeats Location")
            ggtitle("Deletions around repeats (both/one/no breakpoints)")
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_piechart.pdf"),width=5,height=5)

    ## --- count repeat class for homeology length 30bp
    candidates_30bp <- candidates %>% dplyr::filter(min_similarity==80, alignment_length_cutoff %in% c(30), allelic_status_brca1_brca2 != 'other')

    ## for group of allelic_status_brca1_brca2
    del_per_group <- candidates_30bp %>% group_by(allelic_status_brca1_brca2) %>% summarise(n_group=n())

    candidates_30bp_smry <- candidates_30bp %>% group_by(allelic_status_brca1_brca2,repeat_group,query_repClass,subject_repClass) %>% summarise(n=n()) %>% 
        dplyr::mutate(
            query_repClass=ifelse(query_repClass=='','No repeats',query_repClass),
            subject_repClass=ifelse(subject_repClass=='','No repeats',subject_repClass)) %>%
        left_join(del_per_group,by="allelic_status_brca1_brca2") %>%
        dplyr::mutate(prop=n/n_group,repClass=paste(query_repClass,subject_repClass,sep=" <> ")) %>% arrange(desc(prop))
    candidates_30bp_smry$repClass <- factor(candidates_30bp_smry$repClass,levels=as.character(unique(candidates_30bp_smry$repClass)))

    # make stacked barplot
    p <- candidates_30bp_smry %>% ggplot(aes(x=allelic_status_brca1_brca2,y=n,fill=repClass)) + 
        geom_bar(stat="identity",position="fill") + theme_bw() + 
        xlab('') + ylab('Proportion') + 
        ggtitle("Repeat class around deletions with homeology") + scale_fill_hue(name="Repeat Class")
        #theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_repClass_barplot.pdf"),width=7,height=5)

    # make stacked barplot with facet
    p <- candidates_30bp_smry %>% ggplot(aes(fill=allelic_status_brca1_brca2,y=prop,x=repClass)) + 
        geom_bar(stat="identity",position=position_dodge()) + theme_bw() + 
        xlab('') + ylab('Proportion') + 
        ggtitle("Repeat class around deletions with homeology") + scale_fill_viridis_d(name="Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_discrete(drop=F) + facet_wrap(~allelic_status_brca1_brca2)
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_repClass_barplot2.facet.pdf"),width=10,height=5)

    
    ## repClass for all deletions without grouping

    candidates_30bp_all <- candidates_30bp %>% group_by(repeat_group,query_repClass,subject_repClass) %>% summarise(n=n()) %>% 
        dplyr::mutate(
            query_repClass=ifelse(query_repClass=='','No repeats',query_repClass),
            subject_repClass=ifelse(subject_repClass=='','No repeats',subject_repClass)) %>%
        dplyr::mutate(prop=n/nrow(candidates_30bp),repClass=paste(query_repClass,subject_repClass,sep=" <> ")) %>% arrange(desc(prop))
    candidates_30bp_all$repClass <- factor(candidates_30bp_all$repClass,levels=as.character(unique(candidates_30bp_all$repClass)))

    # make pie chart
    p <- candidates_30bp_all %>%
        ggplot(aes(x="",y=prop,fill=repClass)) + geom_bar(stat="identity", width=1) + 
            coord_polar("y", start=0) + 
            theme_void() + scale_fill_hue(name="Repeat Class") +
            ggtitle("Deletions with homeology")
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_repClass_piechart.pdf"),width=7,height=7)

    # make stacked barplot with facet
    p <- candidates_30bp_all %>% ggplot(aes(y=prop,x=repClass)) + 
        geom_bar(stat="identity",position=position_dodge(),fill='royalblue') + theme_bw() + 
        xlab('') + ylab('Proportion') + 
        ggtitle("Repeat class around deletions with homeology") + scale_fill_viridis_d(name="Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_discrete(drop=F)
    ggsave(paste0(outDir,"/repeats_summary.homelen_30_repClass_barplot2.pdf"),width=7,height=5)



}

# the 200 bp flanking sequence around deletions were annotated to repeatmasker
summarize_repeats_for_deletions <- function(deletions_add_repeats,outDir){

    library(tidyr)

    dir.create(outDir, recursive = T, showWarnings = F)

    # group deletions by two breakpoints whether they overlapped with repeats, for example, both ends overlapped with repeats, one end overlapped with repeats, or no end overlapped with repeats
    deletions_add_repeats <- group_deletions_by_repeats(deletions_add_repeats)

    repeats_summary <- deletions_add_repeats %>% 
        #group_by(min_similarity,alignment_length_cutoff,allelic_status_brca1_brca2) %>% 
        summarise(Both=sum(repeat_group=="Both breakpoints"),
                One=sum(repeat_group=="One breakpoint"),
                No=sum(repeat_group=="No repeats"),
                total=nrow(.),
                query_overlap_100pct     =sum(query_rmsk_overlap_width/(query_end-query_start+1)==1),
                query_overlap_above_90pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.9),
                query_overlap_above_80pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.8),
                query_overlap_above_50pct=sum(query_rmsk_overlap_width/(query_end-query_start+1)>0.5),
                subject_overlap_100pct     =sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)==1),
                subject_overlap_above_90pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.9),
                subject_overlap_above_80pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.8),
                subject_overlap_above_50pct=sum(subject_rmsk_overlap_width/(subject_end-subject_start+1)>0.5)
                ) %>%
        dplyr::mutate(pct_Both=Both/total*100,pct_One=One/total*100,pct_No=No/total*100) %>%
        dplyr::mutate(
            pct_query_overlap_100pct     =query_overlap_100pct/(Both+One)*100,
            pct_query_overlap_above_90pct=query_overlap_above_90pct/(Both+One)*100,
            pct_query_overlap_above_80pct=query_overlap_above_80pct/(Both+One)*100,
            pct_query_overlap_above_50pct=query_overlap_above_50pct/(Both+One)*100,
            pct_subject_overlap_100pct     =subject_overlap_100pct/(Both+One)*100,
            pct_subject_overlap_above_90pct=subject_overlap_above_90pct/(Both+One)*100,
            pct_subject_overlap_above_80pct=subject_overlap_above_80pct/(Both+One)*100,
            pct_subject_overlap_above_50pct=subject_overlap_above_50pct/(Both+One)*100
        )
    repeats_summary  %>% fwrite(paste0(outDir,"/repeats_summary.tsv"),sep="\t")

    # make piechart
    p <- repeats_summary %>%
        pivot_longer(cols=c('pct_Both','pct_One','pct_No'),names_to="repeats_location",values_to="pct_repeats") %>%
        dplyr::mutate(repeats_location=case_when(
            repeats_location=="pct_Both" ~ "Both breakpoints",
            repeats_location=="pct_One" ~ "One breakpoint",
            repeats_location=="pct_No" ~ "No repeats",
        )) %>%
        dplyr::mutate(repeats_location=factor(repeats_location,levels=c("Both breakpoints","One breakpoint","No repeats"))) %>%
        ggplot(aes(x="",y=pct_repeats,fill=repeats_location)) + geom_bar(stat="identity", width=1) + 
            coord_polar("y", start=0) + 
            geom_text(aes(label = paste0(round(pct_repeats,0),"%")), size = 6, color = "black",
                position = position_stack(vjust = 0.5)) +
            theme_void() + scale_fill_hue(name="Repeats Location")
            ggtitle("Deletions around repeats (both/one/no breakpoints)")
    ggsave(paste0(outDir,"/repeats_summary.piechart.pdf"),width=5,height=5)

    
    ## make stacked barplot with facet
    deletions_add_repeats_smry <- deletions_add_repeats %>% group_by(repeat_group,query_repClass,subject_repClass) %>% summarise(n=n()) %>% 
        dplyr::mutate(
            query_repClass=ifelse(query_repClass=='' | is.na(query_repClass),'No repeats',query_repClass),
            subject_repClass=ifelse(subject_repClass=='' | is.na(subject_repClass),'No repeats',subject_repClass)) %>%
        dplyr::mutate(prop=n/nrow(deletions_add_repeats),repClass=paste(query_repClass,subject_repClass,sep=" <> ")) %>% arrange(desc(prop))
    deletions_add_repeats_smry$repClass <- factor(deletions_add_repeats_smry$repClass,levels=as.character(unique(deletions_add_repeats_smry$repClass)))
    deletions_add_repeats_smry %>% fwrite(paste0(outDir,"/repeats_summary.repClass.tsv"),sep="\t")

    p <- deletions_add_repeats_smry %>% ggplot(aes(y=prop,x=repClass)) + 
        geom_bar(stat="identity",position=position_dodge(),fill='royalblue') + theme_bw() + 
        xlab('') + ylab('Proportion') + 
        ggtitle("Repeat class around deletions with homeology") + scale_fill_viridis_d(name="Repeat Class") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_x_discrete(drop=F)
    ggsave(paste0(outDir,"/repeats_summary.repClass_barplot.pdf"),width=15,height=5)

}

# summarize repeats in the whole genome
summarize_repeats_in_genome <- function(outDir){

    library(tidyr)
    library(GenomicRanges)

    dir.create(outDir, recursive = T, showWarnings = F)

    # group by repClass and then calculate size of repeats
    rmsk <- fread("data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz") %>%
        dplyr::select(genoName,genoStart,genoEnd,strand,repName,repClass,repFamily)

    repClasses <- unique(rmsk$repClass)

    out <- lapply(repClasses, function(x){

        print(x)

        rmsk_gr <- makeGRangesFromDataFrame(rmsk %>% dplyr::filter(repClass==x) %>% dplyr::mutate(genoName=gsub("chr","",genoName)), 
                seqnames.field = 'genoName', start.field='genoStart', end.field='genoEnd',keep.extra.columns=T)

        red <- reduce(rmsk_gr)
        size <- vapply(width(red), sum, numeric(1))
        sum_size <- sum(size); print(sum_size)
        mean_size <- mean(size)
        return(data.frame(repClass=x,n=length(rmsk_gr),n_collapsed=length(red),sum_size=sum_size,mean_size=mean_size,genome_size=2.45e9,size_prop=sum_size/2.45e9))
    })  

    out <- do.call(rbind,out)

    out %>% fwrite(paste0(outDir,"/repeats_summary_in_genome.tsv"),sep="\t")

    # make barplot
    p <- out %>% dplyr::mutate(repClass=factor(repClass,levels=out %>% arrange(desc(size_prop)) %>% pull(repClass))) %>% 
        ggplot(aes(x=repClass,y=size_prop)) + geom_bar(stat="identity") + 
        theme_bw() + 
        xlab('Repeat Class') + ylab('Proportion in the whole genome') + 
        ggtitle("Repeat class in the whole genome") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0(outDir,"/repeats_summary_in_genome.barplot.pdf"),width=7,height=5)
    
    rm(rmsk)
}