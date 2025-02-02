
library(dplyr)
library(data.table)

prepare_repeatmasker_igv_track <- function(){
    # hg19
    rmsk <- fread("data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz") %>% # /lila/data/riazlab/projects/zhuy1/my_projects/Manisha_SSA/ssa_homeology_20231119/data/RepeatMasker/ucsc_hg19_repeatmasker.txt.gz
        dplyr::rename(chrom=genoName,chromStart=genoStart,chromEnd=genoEnd) %>%
        dplyr::mutate(name=paste0(repClass,"_",repFamily,"_",repName)) %>%
        dplyr::select(chrom,chromStart,chromEnd,name)

    rmsk %>% fwrite("data/RepeatMasker/ucsc_hg19_repeatmasker.bed",sep="\t",col.names=F)

    # hg38
    rmsk <- fread("data/RepeatMasker/ucsc_hg38_repeatmasker.txt.gz") %>% # /lila/data/riazlab/projects/zhuy1/my_projects/Manisha_SSA/ssa_homeology_20231119/data/RepeatMasker/ucsc_hg38_repeatmasker.txt.gz
        dplyr::rename(chrom=genoName,chromStart=genoStart,chromEnd=genoEnd) %>%
        dplyr::mutate(name=paste0(repClass,"_",repFamily,"_",repName)) %>%
        dplyr::select(chrom,chromStart,chromEnd,name)

    rmsk %>% fwrite("data/RepeatMasker/ucsc_hg38_repeatmasker.bed",sep="\t",col.names=F)

}

prepare_homeology_sequence_igv_track <- function(candidates,sample_name=NULL,min_similarity=80,alignment_length_cutoff=30,outDir="igv_tracks"){

    #candidates <- fread("output/serena_SV_deletions_ssa_blast_lr100bp/output_deletions_above_1KB_mostSimilarAlignment/ssa_events_candidates.add_repeats.tsv")
    dir.create(outDir,showWarnings = F,recursive = T)

    min_similarity_thres = min_similarity
    alignment_length_cutoff_thres = alignment_length_cutoff

    if(!is.null(sample_name)){
        selected = candidates %>% 
            dplyr::filter(SAMPLE.TUMOR==sample_name,min_similarity==min_similarity_thres,alignment_length_cutoff==alignment_length_cutoff_thres) %>%
            dplyr::mutate(name=paste0(SAMPLE.TUMOR,"_",CHROM,"_",start_position,"_",end_position)) %>%
            dplyr::select(CHROM,query_start,query_end,CHROM,subject_start,subject_end,name)
    }else{
        selected = candidates %>% 
            dplyr::filter(min_similarity==min_similarity_thres,alignment_length_cutoff==alignment_length_cutoff_thres)

        


        # extract breakpoints
        selected_bedpe = selected %>% 
            dplyr::mutate(name=paste0(SAMPLE.TUMOR,"_",CHROM,"_",start_position,"_",end_position)) %>%
            dplyr::mutate(chrom2=CHROM) %>%
            dplyr::mutate(chrom1=CHROM,start1=start_position-1,end1=start_position,start2=end_position-1,end2=end_position) %>%
            dplyr::select(chrom1,start1,end1,chrom2,start2,end2,name)
        
        selected_bedpe %>% 
            fwrite(paste0(outDir,"/breakpoints_similarity_",min_similarity,"_homeology_length_",alignment_length_cutoff,".bedpe"),sep="\t",col.names=F)

        # extract homeology loci
        selected_homeology_loci_query <- selected %>%
            dplyr::mutate(name=paste0("Query:",SAMPLE.TUMOR,"_",CHROM,"_",query_start,"_",query_end," / ","Subject:",SAMPLE.TUMOR,"_",CHROM,"_",subject_start,"_",subject_end)) %>%
            dplyr::select(CHROM,query_start,query_end,name) %>%
            dplyr::rename(start=query_start,end=query_end)

        selected_homeology_loci_subject <- selected %>%
            dplyr::mutate(name=paste0("Subject:",SAMPLE.TUMOR,"_",CHROM,"_",subject_start,"_",subject_end," / ","Query:",SAMPLE.TUMOR,"_",CHROM,"_",query_start,"_",query_end," / ")) %>%
            dplyr::select(CHROM,subject_start,subject_end,name) %>%
            dplyr::rename(start=subject_start,end=subject_end)

        selected_homeology_loci <- 
            rbind(selected_homeology_loci_query,selected_homeology_loci_subject) 
        
        selected_homeology_loci %>% 
            fwrite(paste0(outDir,"/homeology_loci_similarity_",min_similarity,"_homeology_length_",alignment_length_cutoff,".bed"),sep="\t",col.names=F)


        # convert hg19 to hg38 coordinates
        library(rtracklayer)
        convert_hg19_to_hg38_coordinates <- function(CHROM, start, end){
            brk <- GRanges(seqnames=paste0('chr',CHROM),IRanges(start=start,end=end))
            map.chain <- import.chain("data/hg19ToHg38.over.chain")
            seqlevelsStyle(brk) = "UCSC"
            brk.hg38 <- liftOver(brk, map.chain) # A ‘GRangesList’ object. Each element contains the ranges mapped from the corresponding element in the input (may be one-to-many).
            # brk.hg38 <- unlist(brk.hg38) # directly combine all the ranges may cause problem, because some coordinates may be mapped to multiple locations
            #                                after comparing command line liftOver, the smallest position should be the start, and the largest position should be the end
            brk.hg38.df <- do.call(rbind,lapply(brk.hg38, function(gr){ 
                data.frame(CHROM=seqnames(gr)[1],start_position.hg38=min(start(gr)),end_position.hg38=max(end(gr)))
            }))
            # make sure the order is the same as the original data
            print(all(length(CHROM)==nrow(brk.hg38.df)))
            return(brk.hg38.df)
        }
        # replace hg19 coordinates with hg38 coordinates
        df_hg38 <- convert_hg19_to_hg38_coordinates(selected$CHROM, selected$start_position, selected$end_position)
        selected_hg38 <- selected %>% dplyr::mutate(start_position=df_hg38$start_position.hg38,end_position=df_hg38$end_position.hg38) %>%
            dplyr::select(SAMPLE.TUMOR,CHROM,start_position,end_position)

        df_hg38 <- convert_hg19_to_hg38_coordinates(selected$CHROM, selected$query_start, selected$query_end)
        query_hg38 <- selected %>% dplyr::mutate(query_start=df_hg38$start_position.hg38,query_end=df_hg38$end_position.hg38) %>%
            dplyr::select(SAMPLE.TUMOR,CHROM,query_start,query_end)

        df_hg38 <- convert_hg19_to_hg38_coordinates(selected$CHROM, selected$subject_start, selected$subject_end)
        subject_hg38 <- selected %>% dplyr::mutate(subject_start=df_hg38$start_position.hg38,subject_end=df_hg38$end_position.hg38) %>%
            dplyr::select(SAMPLE.TUMOR,CHROM,subject_start,subject_end)

        selected_hg38 <- cbind(selected_hg38, query_hg38 %>% dplyr::select(query_start,query_end), subject_hg38 %>% dplyr::select(subject_start,subject_end))

        # save breakpoints
        selected_bedpe = selected_hg38 %>% 
            dplyr::mutate(name=paste0(SAMPLE.TUMOR,"_",CHROM,"_",start_position,"_",end_position)) %>%
            dplyr::mutate(chrom2=CHROM) %>%
            dplyr::mutate(chrom1=CHROM,start1=start_position-1,end1=start_position,start2=end_position-1,end2=end_position) %>%
            dplyr::select(chrom1,start1,end1,chrom2,start2,end2,name)
        
        selected_bedpe %>% 
            fwrite(paste0(outDir,"/breakpoints_similarity_",min_similarity,"_homeology_length_",alignment_length_cutoff,".hg38.bedpe"),sep="\t",col.names=F)

        # save homeology loci
        selected_homeology_loci_query <- selected_hg38 %>%
            dplyr::mutate(name=paste0("Query:",SAMPLE.TUMOR,"_",CHROM,"_",query_start,"_",query_end," / ","Subject:",SAMPLE.TUMOR,"_",CHROM,"_",subject_start,"_",subject_end)) %>%
            dplyr::select(CHROM,query_start,query_end,name) %>%
            dplyr::rename(start=query_start,end=query_end)

        selected_homeology_loci_subject <- selected_hg38 %>%
            dplyr::mutate(name=paste0("Subject:",SAMPLE.TUMOR,"_",CHROM,"_",subject_start,"_",subject_end," / ","Query:",SAMPLE.TUMOR,"_",CHROM,"_",query_start,"_",query_end," / ")) %>%
            dplyr::select(CHROM,subject_start,subject_end,name) %>%
            dplyr::rename(start=subject_start,end=subject_end)

        selected_homeology_loci <- 
            rbind(selected_homeology_loci_query,selected_homeology_loci_subject) 
        
        selected_homeology_loci %>% 
            fwrite(paste0(outDir,"/homeology_loci_similarity_",min_similarity,"_homeology_length_",alignment_length_cutoff,".hg38.bed"),sep="\t",col.names=F)


    }
    
}

# to visualize the results in IGV, you need to prepare the tracks
prepare_igv_tracks_marcin <- function(){

    # prepare repeatmasker track
    prepare_repeatmasker_igv_track()

    # prepare breakpoints and homeology loci tracks
    outdir="output/marcin_lr100_blast/output_deletions_above_1KB_mostSimilarAlignment_20240201/igv_tracks"
    candidates <- fread("output/marcin_lr100_blast/output_deletions_above_1KB_mostSimilarAlignment_20240201/ssa_events_candidates.add_repeats.tsv")

    candidates <- candidates %>% filter(min_similarity==80,alignment_length_cutoff==30)

    prepare_homeology_sequence_igv_track(candidates,sample_name=NULL,min_similarity=80,alignment_length_cutoff=30,outDir=outdir)

    # then in IGV, load the tracks
    # 1. load reference genome
    # 2. load repeatmasker track
    # 3. load breakpoints .bedpe tracks, set to "Nested Arcs", "Show Blocks", "Set line thickness to 2"
    # 4. load homeology loci .bed tracks
    # 5. use multi-view to compare the tracks, simply type the two regions in the search box and press enter: 

    # examples: 
    # not clear Query:DO218333_1_11414703_11414796 / Subject:DO218333_1_11417932_11418025
    # Good      Query:1331_11_19245647_19245722 / Subject:1331_11_19260664_19260737
    # Good      Query:1331_13_73240760_73240880 / Subject:1331_13_73246974_73247094
    # not clear Query:1331_20_36475982_36476135 / Subject:1331_20_36482794_36482948
    # Good      Query:1331_3_48143289_48143334 / Subject:1331_3_48145805_48145850
    # Good      Query:1331_7_136762890_136763065 / Subject:1331_7_136769956_136770130

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='LINE',subject_repClass=='LINE')
    # 10 deletions
    # 1331 chr7:136447736-136454780 chr7:136446736-136448736 chr7:136453780-136455780
    # hg38: chr7:136762989-136770033

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='LTR',subject_repClass=='LTR')
    # 6 deletions

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='SINE',subject_repClass=='SINE')
    # 184 deletions

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='Simple_repeat'|subject_repClass=='Simple_repeat')) %>% 
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass)

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='Low_complexity'|subject_repClass=='Low_complexity')) %>%
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass)

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='DNA'|subject_repClass=='DNA')) %>%
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass) %>% head
    
}


# to visualize the results in IGV, you need to prepare the tracks
prepare_igv_tracks_serena_breast <- function(){

    # prepare repeatmasker track
    prepare_repeatmasker_igv_track()

    outdir="output/serena_SV_deletions_ssa_blast_lr100bp/output_deletions_above_1KB_mostSimilarAlignment_updatedControls/igv_tracks"

    # prepare breakpoints and homeology loci tracks
    candidates <- fread("output/serena_SV_deletions_ssa_blast_lr100bp/output_deletions_above_1KB_mostSimilarAlignment_updatedControls/ssa_events_candidates.add_repeats.tsv")

    candidates <- candidates %>% filter(min_similarity==80,alignment_length_cutoff==30)

    prepare_homeology_sequence_igv_track(candidates,sample_name=NULL,min_similarity=80,alignment_length_cutoff=30,outDir="igv_tracks")

    # then in IGV, load the tracks
    # 1. load reference genome
    # 2. load repeatmasker track
    # 3. load breakpoints .bedpe tracks, set to "Nested Arcs", "Show Blocks", "Set line thickness to 2"
    # 4. load homeology loci .bed tracks
    # 5. use multi-view to compare the tracks, simply type the two regions in the search box and press enter: 
    # chr2:101924613-101926613 chr2:101983606-101985606

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='LINE',subject_repClass=='LINE')
    # 32 deletions
    # chr7:119850088-119852088 chr7:119851626-119853626

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='LTR',subject_repClass=='LTR')
    # 11 deletions
    # chr5:60500376-60502376 chr5:60502298-60504298

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,query_repClass=='SINE',subject_repClass=='SINE')
    # 400 deletions
    # chr22:40917683-40919683 chr22:40924475-40926475

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='Simple_repeat'|subject_repClass=='Simple_repeat')) %>% 
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass)
    # 2:211163030-211163091 2:211244076-211244109

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='Low_complexity'|subject_repClass=='Low_complexity')) %>%
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass)
    # 5:95967729-95967752 5:101658561-101658595

    candidates %>% filter(min_similarity==80,alignment_length_cutoff==30,(query_repClass=='DNA'|subject_repClass=='DNA')) %>%
        select(CHROM,query_rmsk_start,query_rmsk_end,subject_rmsk_start,subject_rmsk_end,query_repClass,subject_repClass) %>% head
    # 1:191745241-191746542 1:191764672-191765212
    
}