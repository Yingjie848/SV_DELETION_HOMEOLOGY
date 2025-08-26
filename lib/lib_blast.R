
library(Biostrings)

b37 <- readDNAStringSet("~/my_databases/Pier_reference_b37/b37.cleanName.fasta")

# This is the final version of the function to extract flanking sequences of breakpoints
extract_sequence_lr100 <- function(sampleid,chr,start,end,fout_left,fout_right){

    # left 100 bp, including the first breakpoint, right 100 bp, not including the first breakpoint
    lseq=gsub("^ ","",toString(subseq(b37[chr],start=start-99,end=start+100)))
    # left 100 bp not including the second breakpoint, right 100 bp, including the second breakpoint
    rseq=gsub("^ ","",toString(subseq(b37[chr],start=end-100,end=end+99)))

    eventid <- paste0(sampleid,"/",chr,':',start,"-",end)

    lseqname=paste0(">",eventid,",",start-99,"-",start+100)
    rseqname=paste0(">",eventid,",",end-100  ,"-",end+99)

    cat(lseqname,"\n",lseq,"\n",file=fout_left,append=TRUE,sep="")
    cat(rseqname,"\n",rseq,"\n",file=fout_right,append=TRUE,sep="")

}


# extract flacking sequences of breakpoints just for later checking
extract_flanking_sequence_serena_data <- function(d.out,blast_dir){

    #library(Biostrings)
    #b37 = readDNAStringSet("~/my_databases/Pier_reference_b37/b37.cleanName.fasta")

    # example: PD8619, chr11, 67851596, 67855463, 3868, CATGCCTGTAATCCCAGCTACTCAGG
    sampleid='PD8619'; chr=11; start=67851596; end=67855463; fout_left=paste0(blast_dir,"/lseq.fasta"); fout_right=paste0(blast_dir,"/rseq.fasta")
    extract_sequence_lr100(sampleid,chr,start,end,fout_left,fout_right)

    # extract all events just for later checking
    fout_left =paste0(blast_dir,"/all_lseq.fasta")
    fout_right=paste0(blast_dir,"/all_rseq.fasta")
    if(file.exists(fout_left)){
        file.remove(fout_left)
    }
    if(file.exists(fout_right)){
        file.remove(fout_right)
    }

    out <- lapply(1:nrow(d.out),function(i){
        sampleid <- d.out[i,]$SAMPLE.TUMOR
        chr      <- d.out[i,]$CHROM
        start    <- d.out[i,]$start_position
        end      <- d.out[i,]$end_position
        #cat(sampleid,chr,start,end,"\n")

        extract_sequence_lr100(sampleid,chr,start,end,fout_left,fout_right)
    })

}

# extract flacking sequences of breakpoints just for later checking
extract_flanking_sequence_marcin_data <- function(d.out,blast_dir){

    #library(Biostrings)
    #b37 = readDNAStringSet("~/my_databases/Pier_reference_b37/b37.cleanName.fasta")

    # extract all events
    fout_left =paste0(blast_dir,"/all_lseq.fasta")
    fout_right=paste0(blast_dir,"/all_rseq.fasta")
    if(file.exists(fout_left)){
        file.remove(fout_left)
    }
    if(file.exists(fout_right)){
        file.remove(fout_right)
    }


    out <- lapply(1:nrow(d.out),function(i){
        sampleid <- d.out[i,]$SAMPLE.TUMOR
        chr      <- d.out[i,]$CHROM
        start    <- d.out[i,]$start_position
        end      <- d.out[i,]$end_position
        #cat(sampleid,chr,start,end,"\n")

        extract_sequence_lr100(sampleid,chr,start,end,fout_left,fout_right)
    })

}




# This is the final version of the function to run pairwise blast for flanking +/- 100 bp sequences
# It calls blastn on server and submit job to clusters via sub.sh
run_pairwise_blast_lr100 <- function(deletion_tbl,bltbl_file,tmp_dir="tmp"){

    samples <- unique(deletion_tbl$SAMPLE.TUMOR)

    unlink(tmp_dir,recursive=TRUE)
    file.remove(paste0(tmp_dir,".run.sh"))
    file.remove(paste0(tmp_dir,".run.err"))
    file.remove(paste0(tmp_dir,".run.out"))

    out <- lapply(1:nrow(deletion_tbl),function(i){

        query <- deletion_tbl[i,]$SAMPLE.TUMOR; cat(i," ",query,"\n")

        tmp_subdir <- paste0(tmp_dir,"/",query); dir.create(tmp_subdir,recursive=TRUE)

        sampleid <- gsub("/.*","",deletion_tbl[i,]$SAMPLE.TUMOR)
        chr      <- gsub(".*/","",gsub(":.*","",deletion_tbl[i,]$CHROM))
        start    <- as.numeric(gsub("-.*","",gsub(".*:","",deletion_tbl[i,]$start_position)))
        end      <- as.numeric(gsub(".*-","",gsub(",.*","",deletion_tbl[i,]$end_position)))
        #cat(sampleid,chr,start,end,"\n")

        fout_left =paste0(tmp_subdir,"/del_",i,"_lseq.fasta")
        fout_right=paste0(tmp_subdir,"/del_",i,"_rseq.fasta")
        extract_sequence_lr100(sampleid,chr,start,end,fout_left,fout_right)

        # build database
        cmd1 <- paste0("makeblastdb -in ",fout_right," -dbtype nucl; "); #system(cmd)
        #cat(cmd,file=paste0(tmp_subdir,".run.sh"),append=TRUE)
        #cat(cmd,file=paste0(tmp_dir,".run.sh"),append=TRUE)

        # run blastn
        bltbl <- paste0(tmp_subdir,"/del_",i,".blast.tsv")
        cmd2 <- paste0("blastn -query ",fout_left," -db ",fout_right," -out ",bltbl," -word_size 4 -evalue 1000 -outfmt 6 -dust no -soft_masking false \n"); # system(cmd)
        #cat(cmd2,file=paste0(tmp_subdir,".run.sh"),append=TRUE)
        cat(cmd1,cmd2,file=paste0(tmp_dir,".run.sh"),append=TRUE)

    })

    # submit jobs by patients
    #out <- sapply(samples,function(sample){
    #    err_file=paste0(tmp_dir,"/",sample,".run.err")
    #    out_file=paste0(tmp_dir,"/",sample,".run.out")
    #    cmd <- paste0("bsub -W 2:00 -R 'rusage[mem=1]' -e ",err_file," -o ",out_file," sh ",tmp_dir,"/",sample,".run.sh"); system(cmd)
    #})

    err_file=paste0(tmp_dir,".run.err")
    out_file=paste0(tmp_dir,".run.out")
    cmd <- paste0("lib/sub.sh 12:00 1 1 ",err_file," ",out_file,' "sh ',tmp_dir,'.run.sh"'); system(cmd)

    # wait until all jobs finished
    cmd <- paste0("find ",tmp_dir," | awk '$1~/blast.tsv/' | awk '{print \"cat \"$1}' | bash > ",bltbl_file); system(cmd)

}


# Process blast output in blast_dir/blast_output.tsv
# add columns: SAMPLE.TUMOR, CHROM, start_position, end_position, q_strand, s_strand, query_key, subject_key
# save back to blast_dir/blast_output.tsv
process_blast_output <- function(blast_dir){

    blast_output <- fread(paste0(blast_dir, "/blast_output.tsv"))
    names(blast_output) <- c("query", "subject", "identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

    blast_output <- blast_output %>%
        dplyr::mutate(query_key = gsub(",.*", "", query), subject_key = gsub(",.*", "", subject)) %>%
        dplyr::filter(query_key == subject_key) %>%
        dplyr::mutate(
            q_strand = ifelse(q_start <= q_end, "+", "-"),
            s_strand = ifelse(s_start <= s_end, "+", "-")
        )

    split_query_id <- lapply(blast_output$query, function(x) strsplit(x, "/|:|-|,"))

    blast_output <- blast_output %>%
        dplyr::mutate(
            SAMPLE.TUMOR = sapply(1:length(split_query_id), function(i) split_query_id[[i]][[1]][1]),
            CHROM = sapply(1:length(split_query_id), function(i) split_query_id[[i]][[1]][2]),
            start_position = sapply(1:length(split_query_id), function(i) as.numeric(split_query_id[[i]][[1]][3])),
            end_position = sapply(1:length(split_query_id), function(i) as.numeric(split_query_id[[i]][[1]][4]))
        )

    # write out formatted blast results
    blast_output %>% fwrite(paste0(blast_dir,"/blast_output.formatted.tsv"), sep = "\t")

    blast_output_same_strand <- blast_output %>% dplyr::filter(q_strand=='+',s_strand=='+')
    blast_output_same_strand %>% fwrite(paste0(blast_dir,"/blast_output.same_direction.tsv"), sep = "\t")

}

# Process blast output in blast_dir/blast_output.tsv
# add columns: SAMPLE.TUMOR, CHROM, start_position, end_position, q_strand, s_strand, query_key, subject_key
# save back to blast_dir/blast_output.tsv
process_blast_output_serena_data <- function(blast_dir){

    blast_output <- fread(paste0(blast_dir,"/blast_output.tsv"))
    names(blast_output) <- c('query','subject','identity','alignment_length','mismatches','gap_opens','q_start','q_end','s_start','s_end','evalue','bit_score')

    blast_output <- blast_output %>% 
                    dplyr::mutate(query_key=gsub(",.*","",query), subject_key=gsub(",.*","",subject)) %>% 
                    dplyr::filter(query_key==subject_key) %>%
                    dplyr::mutate(q_strand=ifelse(q_start<=q_end,'+','-'),
                                s_strand=ifelse(s_start<=s_end,'+','-'))

    split_query_id <- lapply(blast_output$query,function(x) strsplit(x,"/|:|-|,"))

    blast_output <- blast_output %>%
                    dplyr::mutate(SAMPLE.TUMOR = sapply(1:length(split_query_id),function(i) split_query_id[[i]][[1]][1]),
                                CHROM        = sapply(1:length(split_query_id),function(i) split_query_id[[i]][[1]][2]),
                                start_position = sapply(1:length(split_query_id),function(i) as.numeric(split_query_id[[i]][[1]][3])),
                                end_position   = sapply(1:length(split_query_id),function(i) as.numeric(split_query_id[[i]][[1]][4])))

    # write out formatted blast results
    blast_output %>% fwrite(paste0(blast_dir,"/blast_output.formatted.tsv"), sep = "\t")

    blast_output_same_strand <- blast_output %>% dplyr::filter(q_strand=='+',s_strand=='+')
    blast_output_same_strand %>% fwrite(paste0(blast_dir,"/blast_output.same_direction.tsv"), sep = "\t")

}

# Process blast output in blast_dir/blast_output.tsv
# add columns: SAMPLE.TUMOR, CHROM, start_position, end_position, q_strand, s_strand, query_key, subject_key
# save back to blast_dir/blast_output.tsv
process_blast_output_marcin_data <- function(blast_dir){
    
    blast_output <- fread(paste0(blast_dir,"/blast_output.tsv"))
    names(blast_output) <- c('query','subject','identity','alignment_length','mismatches','gap_opens','q_start','q_end','s_start','s_end','evalue','bit_score')

    blast_output <- blast_output %>% 
                    dplyr::mutate(query_key=gsub(",.*","",query), subject_key=gsub(",.*","",subject)) %>% 
                    dplyr::filter(query_key==subject_key) %>%
                    dplyr::mutate(q_strand=ifelse(q_start<=q_end,'+','-'),
                                s_strand=ifelse(s_start<=s_end,'+','-')) 

    # write out formatted blast results
    blast_output %>% fwrite(paste0(blast_dir,"/blast_output.formatted.tsv"), sep = "\t")


    ######################################################################################################
    # extract alignments aligned on the same strand

    blast_output_same_strand <- blast_output %>% dplyr::filter(q_strand=='+',s_strand=='+')

    # to avoid conflict of '-' in sample name, replace it with '__', then replace it back

    split_query_id <- lapply(blast_output_same_strand$query,function(x) strsplit(x,"/"))

    SAMPLE.TUMOR_   = sapply(1:length(split_query_id),function(i) split_query_id[[i]][[1]][1])

    split_query_id <- lapply(sapply(1:length(split_query_id),function(i) split_query_id[[i]][[1]][2]),function(x) strsplit(x,":|-|,"))

    CHROM_          = sapply(1:length(split_query_id),function(i) split_query_id[[i]][[1]][1])
    start_position_ = sapply(1:length(split_query_id),function(i) as.numeric(split_query_id[[i]][[1]][2]))
    end_position_   = sapply(1:length(split_query_id),function(i) as.numeric(split_query_id[[i]][[1]][3]))

    blast_output_same_strand <- blast_output_same_strand %>%  
                                dplyr::mutate(  SAMPLE.TUMOR   = SAMPLE.TUMOR_,
                                                CHROM          = CHROM_,
                                                start_position = start_position_,
                                                end_position   = end_position_) 
    blast_output_same_strand %>% fwrite(paste0(blast_dir,"/blast_output.same_direction.tsv"), sep = "\t")

}


