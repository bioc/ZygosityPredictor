#' @keywords internal
catt <- function(printLog, level, text){
  if(printLog==TRUE){
    message(rep("  ", level), text)
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr tibble
prepare_raw_bam_file <- function(bamDna, chr1, chr2, pos1, pos2, 
                                 verbose=FALSE){
  qname.first <- . <- NULL
  ## importFrom dplyr tibble filter
  ref_pos1 <- as.numeric(pos1)
  ref_pos2 <- as.numeric(pos2)
  ref_chr1 <- as.character(chr1)
  ref_chr2 <- as.character(chr2)
  if(pos1>pos2){
    ref_pos1 <- as.numeric(pos2)
    ref_pos2 <- as.numeric(pos1)
    ref_chr1 <- as.character(chr2)
    ref_chr2 <- as.character(chr1)
  }
  ref_gr1 <- GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1)
  ref_gr2 <- GRanges(seqnames = ref_chr2, 
                                   ranges = ref_pos2)
  ## now load all reads/read-pairs that cover the position of the first variant
  #vm("loading reads", verbose)
  all_covering_read_pairs <- readGAlignmentPairs(
    bamDna,
    param=ScanBamParam(
      which=GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1),
      what=c("qname","seq", "cigar", "mapq", "qual")
    )) 
  if(length(all_covering_read_pairs)==0){
    return(tibble())
  } else {
    ## combine all ranges and check for ref_pos2
    all_reads <- c(
      GenomicAlignments::first(all_covering_read_pairs) %>%
        GRanges(),
      GenomicAlignments::last(all_covering_read_pairs) %>%
        GRanges()
    )
    
    shared_read_pairs <- all_reads %>%
      subsetByOverlaps(.,ref_gr2) %>%
      elementMetadata(.) %>%
      .[["qname"]] %>%
      unique()
    
    filtered_reads <- all_reads[
      which(all_reads$qname %in% shared_read_pairs)
    ]
    return(filtered_reads)
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom purrr compact
#' @importFrom dplyr summarize pull
insert_missing_cnv_regions <- function(somCna, geneModel, sex, ploidy){
  seqnames <- start <- end <- tcn <- . <- NULL
  if(sex=="male"){
    allowed_chr <- allowed_inputs("chrom_names")
  } else {
    allowed_chr <- allowed_inputs("chrom_names") %>% 
      .[which(!str_detect(., "Y"))]
  }
  tbl_refgen <-
    as_tibble(geneModel) %>%
    group_by(seqnames) %>%
    summarize(mn=min(start),
              mx=max(end)) %>%
    filter(seqnames %in% allowed_chr)
  tbl_cnv <- 
    as_tibble(somCna) %>%
    mutate(tcn=as.numeric(tcn))
  new_somCna <- lapply(as.character(tbl_refgen$seqnames), function(CHR){
    all_chr <- tbl_cnv %>%
      filter(seqnames==CHR) %>%
      arrange(start)
    mn <- tbl_refgen %>% 
      filter(seqnames==CHR) %>% pull(mn)
    mx <- tbl_refgen %>% 
      filter(seqnames==CHR) %>% pull(mx)    
    
    if(nrow(all_chr)==0){
      tbl_refgen %>% 
        filter(seqnames==CHR) %>%
        select(seqnames, start=mn, end=mx) %>% 
        mutate(tcn=ploidy, cna_type="no CNA type annotated") %>%
        GRanges() %>%
        return()
    } else {
      if(all_chr[1,]$start<mn){mn <- all_chr[1,]$start}
      if(all_chr[nrow(all_chr),]$end>mx){mx <- all_chr[nrow(all_chr),]$end}
      poslist <- 
        apply(all_chr, 1, function(ROW){
          return(c(start=as.numeric(ROW[["start"]]), 
                   end=as.numeric(ROW[["end"]])))
        }, simplify=FALSE) %>%
        append(list(c(start=NA, end=mn)),.) %>%
        append(.,list(c(start=mx, end=NA)))
      
      to_insert <- lapply(c(2:length(poslist)), function(I){
        diff <- poslist[[I]][["start"]]-poslist[[I-1]][["end"]]
        if(diff>1){
          return(c(start=poslist[[I-1]][["end"]]+1, 
                   end=poslist[[I-1]][["end"]]+diff-1))
        }
      })  %>%
        compact()
      
      if(length(to_insert)==0){
        GRanges(all_chr) %>% return()
      } else {
        return(c(GRanges(all_chr),  bind_rows(to_insert) %>%
                   mutate(chr=CHR, tcn=ploidy, 
                          cna_type="no CNA type annotated") %>%
                   GRanges()))        
      }
    }
  })
  comb_somCna <- Reduce(function(x,y)c(x,y),new_somCna)
  return(comb_somCna)
}

calc_left_wt_copies <- function(mtcn, status, aff_copies1, aff_copies2){
  left_wt_copies <- ifelse(status=="diff",
                           mtcn-sum(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)),
                           mtcn-max(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)))
  return(as.numeric(left_wt_copies))
}

pre_scoring <- function(tcn1, tcn2, status, aff_copies1, aff_copies2){
  mtcn <- min(as.numeric(tcn1), as.numeric(tcn2))
  left_wt_copies <- calc_left_wt_copies(mtcn, status, aff_copies1, aff_copies2)
  pre_score <- ifelse(status=="diff"&left_wt_copies<0.5,
                      2, 1)
  return(list(mtcn=mtcn,
              left_wt_copies=left_wt_copies,
              pre_score=pre_score))
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr as_tibble group_by mutate tibble tally 
#' @importFrom purrr set_names
classify_combination <- function(classified_reads, purity, eval_full, printLog,
                                 tcn1, tcn2, aff_copies1, 
                                 aff_copies2){
  result <- . <- fac <- NULL
  if(nrow(classified_reads)==0){
    return(rep(0, 7) %>% t() %>% as.data.frame() %>% 
             set_names(c("both", "none", "mut1", "mut2", "dev_var", 
                         "no_overlap", "none_raw")) %>% 
             mutate(status="null"))
  }
  number <- classified_reads %>%
    mutate(fac=factor(result, levels = c('both', 'mut1', 'mut2', 'none', 
                                         "dev_var", 'read_in_read', 
                                         'spanned_out'))) %>%
    group_by(fac, .drop = FALSE) %>%
    tally() %>% 
    column_to_rownames(var='fac') %>%
    t() %>%
    as_tibble()
  both <- number[['both']]
  mut1 <- number[['mut1']]
  mut2 <- number[['mut2']]
  none_raw <- number[['none']]
  # now the purity is calculated out... actually all "none" detetions  count 
  # for mono-allelic inactivation.... but to get a more reliable result it 
  # would be necessary to detect reads which have both mutations
  # therefore we check if any "both" detections were found and if not, we 
  # assume that the "none" detections are artefacts or non tumor cell tissue   
  none_adj_im <- none_raw-sum(number[1,])*(1-as.numeric(purity))
  none_adj <- ifelse(none_adj_im<0,
                     0,
                     round(none_adj_im))
  clear_sum <- sum(none_adj, both, mut1, mut2)
  prob_diff <- sum(mut1, mut2)/clear_sum
  prob_same <- sum(none_adj, both)/clear_sum
  ## here the main decision is made
  if(clear_sum==0){
    status <- "null"
    status_def_info <- "non-overlapping reads only -> null"
  } else if(sum(mut1, mut2, both)==0){
    ## if no evidence for both muts at the same read or at different reads can
    ## be found, we assign status "null"
    status <- "null"
    status_def_info <- 
      "Only reads without any variant: No evidence for same or diff -> null"
#  } else if(eval_full==FALSE&none_raw==0&(mut1==0|mut2==0)){
    ## if we are at extended SNP phasing, indicated by the eval_full parameter
    ## we have to check if the variant we use to phase is a homozygous germline
    ## variant. In this case, we cannot use it for phasing, as it will be 
    ## detected at every read. To filter out such variants, we check if there 
    ## are any reads without the variant. If we cannot find any
    ## read in SNP phasing not having the variant, we assign "null" status as
    ## we assume the variant to phase is homozygous germline. this method is
    ## conservative.
#    status <- "null"
#    status_def_info <- paste(
#      "one SNP of combination seems to be germline homozygous:",
#      "Not usable for phasing -> null")
  } else if(both==0){
    status <- "diff"
    status_def_info <- ifelse(mut1==0|mut2==0,
      paste("no reads with both variants and reads with only one -> diff;",
      "but one variant without any detection!"),
      "variants found only on different reads -> diff"
    )
  } else if(prob_same>prob_diff){
    status <- "same"
    status_def_info <- ifelse(prob_diff>0,
      paste("Reads detected carrying both variants -> same;",
            "but also reads with only one!", "prob same:", prob_same, 
            "prob diff:", prob_diff),
      "variants detected on the same reads -> same"
      )
  } else if((mut1<0.1*sum(mut1, mut2)|mut2<0.1*sum(mut1, mut2))&
            both>0.1*sum(mut1, mut2)){
    status <- "same"
    status_def_info <- 
      paste("reads found with both variants and also with only one",
            "number of detections of one variant smaller than 10 % of",
            "overall evidence for diff and evidence for same higher -> same")
  } else {
    status <- "diff"
    status_def_info <- "variants found only on different reads -> diff"
  }
  ntabs <- ifelse(eval_full==TRUE, 3, 7)
  catt(printLog, ntabs, status_def_info)
  status_table <- tibble(
    both=both,
    none=none_adj,
    mut1=mut1,  
    mut2=mut2,
    dev_var=number[['dev_var']],
    no_overlap=sum(as.numeric(number[['read_in_read']]), 
                   as.numeric(number[['spanned_out']])),
    status=status,
    none_raw=none_raw
    )
  if(eval_full==TRUE){
    psc <- pre_scoring(tcn1, tcn2, status, aff_copies1, aff_copies2)
    if(status=="diff"&psc$pre_score==1&round(psc$mtcn)<=2){
      catt(printLog, 3, 
           "muts at different copies and tcn<=2, but left wt-copies >= 0.5!")
    }
    status_table <- status_table %>%
      mutate(
        score=psc$pre_score,
        wt_cp=round(psc$left_wt_copies,2),
        tcn=round(psc$mtcn,2),
        info=paste(status_def_info, "-> left wt copies", 
                   round(psc$left_wt_copies,2))
      )
  } 
  return(status_table)
}  
#' @keywords internal
#' @importFrom stringr %>% str_detect str_match_all str_sub
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate left_join bind_rows mutate tibble filter
make_read_overview <- function(read_start, seq, cigar, qual){
  . <- element <- width <- end <- NULL
  ## first probably replaceable by cigarRangesAlongQuerySpace(cigar)
  raw_cigs <- str_match_all(cigar,
                            '\\d*[:upper:]') %>%
    unlist() %>%
    tibble(element=.) %>%
    mutate(
      width=str_sub(element, start = 1, end=nchar(element)-1) %>%
        as.numeric(),
      type=str_sub(element, start = nchar(element), end=nchar(element))
    ) %>%
    filter(width!=0) %>%
    rownames_to_column("id")
  raw_maps <-
    raw_cigs[which(!str_detect(raw_cigs$element, "D|N")),] %>%
    mutate(to_extract=cumsum(width)) %>%
    apply(.,1,function(cig){
      if(cig[["id"]]==1){
        curr_seq <- str_sub(seq, start=1, end=cig[["width"]])
        curr_qual <- str_sub(qual, start=1, end=cig[["width"]])
      } else {
        curr_seq <- str_sub(seq,
                            start=as.numeric(cig[["to_extract"]])-
                              as.numeric(cig[["width"]])+1,
                            end=cig[["to_extract"]])
        curr_qual <- str_sub(qual,
                            start=as.numeric(cig[["to_extract"]])-
                              as.numeric(cig[["width"]])+1,
                            end=cig[["to_extract"]])
      }
      return(list(id=cig[["id"]],
                  seq=curr_seq,
                  qual=curr_qual))
    }) %>% bind_rows()
  comb <- left_join(raw_cigs, raw_maps, by="id")
  rel <-
    comb[which(comb$type %in% c("M", "D", "N")),] %>%
    mutate(end=cumsum(width),
           start=end-width+1,
           map_start=read_start+end-width,
           map_end=read_start+end-1) %>%
    select(id, map_start, map_end, start, end)
  full <- left_join(comb, rel, by="id")
  return(full)
}
make_read_overview_new <- function(read_start, seq, cigar, qual){
    ## parse cigar string according to query
  cigq <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar,
                                                        with.ops = T)%>%
    unlist() %>% as.data.frame() %>%
    set_names(.,nm=paste0(names(.), "_q"))
  ## parse cigar string according to reference
  cigr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar,
                                                        with.ops = F)%>%
    unlist() %>% as.data.frame() %>%
    set_names(.,nm=paste0(names(.), "_r"))
  ## combine them and extract relevant bases and quality from sequence
  raw_cigs_new <- bind_cols(cigq, cigr) %>%
    rowwise() %>%
    mutate(seq=str_sub(seq, start=start_q, end=end_q) %>% na_if(""),
           qual=str_sub(qual, start=start_q, end=end_q) %>% na_if(""),
           width=max(width_q, width_r),
           map_start=read_start+end_r-width,
           map_end=read_start+end_r-1) %>%
    select(width, type=names_q, seq, qual, map_start, map_end, start=start_r, 
           end=end_r)
  return(raw_cigs_new)
}

extract_subseq <- function(ref_pos, element_mut, parsed_read, 
                           length_indel, string){
  # print(parsed_read)
  # print(str_sub(element_mut[[string]], start=ref_pos-element_mut$map_start+1, 
  #               end=-1))
  # print(      paste(parsed_read[seq(as.numeric(element_mut$id)+1, 
  #                                   nrow(parsed_read), 1)][[string]], collapse=""))
  pos_in_seq <- ref_pos-element_mut$map_start+1
  subseq_in_element <- str_sub(element_mut[[string]], start=pos_in_seq, 
                               end=-1)
  subseq_in_following_elements <- parsed_read[which(as.numeric(parsed_read$id)>as.numeric(element_mut$id)),] %>%
    pull(all_of(string)) %>%
    na.omit() %>%
    paste(collapse="")
  
  full_remaining <- paste0(subseq_in_element, subseq_in_following_elements)
 # print(parsed_read)
#  print(subseq_in_element)
#  print(subseq_in_following_elements)
#  print(full_remaining)
 # print(length_indel)
  
  return(str_sub(full_remaining, start=1, end=length_indel))
}


extract_snv <- function(ref_pos, element_mut){
  ## get cigar element in which reference position in located
  pos_in_seq <- ref_pos-element_mut$map_start+1
    base <- 
      str_sub(element_mut$seq,
              start=pos_in_seq,
              end=pos_in_seq) 
    qual <- 
      str_sub(element_mut$qual,
              start=pos_in_seq,
              end=pos_in_seq)  
  return(tibble(base=base, qual=qual, info="mapped", len_indel=NA,
                exp_indel=NA))
}

extract_insertion <- function(ref_pos, parsed_read, element_mut, length_ins){
  ## insertion are always indicated by the I type. The reference position tells
  ## the position in front of the I segment, which means we are looking for
  ## the subsequent one type
  expected_indel_detected <- FALSE
  len <- NA
  if(as.numeric(element_mut$id)==nrow(parsed_read)){
    if(ref_pos==element_mut$map_end){
      info <- 
        "no insertion detected: ref pos is last mapped base of read" 
    } else {
      info <- 
        "no insertion detected: ref pos mapped in last element of parsed read"      
    }
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  } else if("I" %in% parsed_read$type){
    ## insertion detected... check if it is the matching one
    if(ref_pos==element_mut$map_end){
      ## if insertion is at refrence position, the subsequent position is an a
      ## new elemnt labelled as I
      next_element <- parsed_read[as.numeric(element_mut$id)+1,]
      if(next_element$type=="I"){
        ## insertion detected at position.. extract inserted sequence
        base <- paste0(
          str_sub(element_mut$seq, start=-1, end=-1),
          next_element$seq
        )
        qual <- paste0(
          str_sub(element_mut$qual, start=-1, end=-1),
          next_element$qual
        )
        info <- "insertion detected"
        len <- next_element$width
        expected_indel_detected <- TRUE
      } else {
        info <- paste(next_element$type, "detected - wrong class")
        base <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, "seq")
        qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, "qual")
      }
    } else {
      info <- "no insertion detected: ref pos not mapped to end of element"
      base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
      qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
    }
  } else {
    ## no inserttion detected.. extract base at position to get base quality
    info <- "no insertion detected: no I in cigar"
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  }
  return(tibble(base=base, qual=qual, info=info, len_indel=len, 
                exp_indel=expected_indel_detected))
}

extract_deletion <- function(ref_pos, parsed_read, element_mut, length_del){
  expected_indel_detected <- FALSE
  len <- NA
  if(as.numeric(element_mut$id)==nrow(parsed_read)){
    if(ref_pos==element_mut$map_end){
      info <- 
        "no deletion detected: ref pos is last mapped base of read" 
    } else {
      info <- 
        "no delection detected: ref pos mapped in last element of parsed read"      
    }
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  } else if("D" %in% parsed_read$type){
    ## insertion detected... check if it is the matching one
    if(ref_pos==element_mut$map_end){
      ## if insertion is at refrence position, the subsequent position is an a
      ## new elemnt labelled as D
      next_element <- parsed_read[as.numeric(element_mut$id)+1,]
      if(next_element$type=="D"){
        expected_indel_detected <- TRUE
        info <- "deletion detected"
        len <- next_element$width
      } else {
        info <- "no deletion detected: subsequent element not I"
      }
    } else {
      info <- "no insertion detected: ref pos not mapped to end of element"
    }      
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, "qual")
  } else {
    ## no inserttion detected.. extract base at position to get base quality
    info <- "no insertion detected: no I in cigar"
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, "qual")
  }
  return(tibble(base=base, qual=qual, info=info, len_indel=len,
                exp_indel=expected_indel_detected))
}

cigar_element <- function(parsed_read, ref_pos){
  parsed_read %>% rowwise() %>%
    filter(between(ref_pos, map_start, map_end)) %>%
    .[1,] %>% return()
}

extract_base_at_refpos <- function(parsed_read, ref_pos, class, ref_alt, ref_ref){
  ## extract elemnt containing reference position from parsed read
  element_mut <- cigar_element(parsed_read, ref_pos)
  ## if position is inside N or D type, the read does not cover the position
  ## N means for example skipped exon in RNA
  ## D wouldmean the position is deleted
  ## real deletions are still mapped in the M part
  if(element_mut$type %in% c("N", "D")){
    base_info <- tibble(base=NA, qual=NA, info="position skipped/deleted", 
                        len_indel=NA, exp_indel=NA)
  } else if(class=="snv"){
    base_info <- extract_snv(ref_pos, element_mut)
  } else if(class=="ins"){
    base_info <- extract_insertion(ref_pos, parsed_read, element_mut, 
                                   nchar(ref_alt)-nchar(ref_ref))
  } else if(class=="del"){
    base_info <- extract_deletion(ref_pos, parsed_read, element_mut, 
                                  1+abs(nchar(ref_ref)-nchar(ref_alt)))
  } else {
    stop("class of variant needs to be provided")
  }
  return(base_info %>% mutate(class=class, mapq=element_mut$mapq))
}

#' @keywords internal
#' @importFrom stringr %>% str_replace str_sub
#' @importFrom dplyr between case_when rowwise filter
check_mut_presence <- function(read_structure, ref_pos,
                                  ref_alt, ref_ref, ref_class, verbose=FALSE){
  map_start <- map_end <- . <- NULL
  ## check if read is spanned out:
  element_mut <- read_structure %>% rowwise() %>%
      filter(between(ref_pos, map_start, map_end)) %>%
      .[1,]
  if(element_mut$type %in% c("N", "D")){
    mut_at_read <- -10
  } else {
    if(ref_class=="snv"){
      pos_in_seq <- ref_pos-element_mut$map_start+1
      base <- 
        str_sub(element_mut$seq,
                start=pos_in_seq,
                end=pos_in_seq) 
      qual <- 
        str_sub(element_mut$qual,
                start=pos_in_seq,
                end=pos_in_seq)
      #vm(paste("base:",base), verbose)
      #vm(paste("alt:",ref_alt), verbose)
      #vm(paste("ref:",ref_ref), verbose)
      mut_at_read <- case_when(
        base==ref_alt ~ 1,
        base==ref_ref ~ 0,
        TRUE ~ -2
      ) 
      #vm(paste("numcat:",mut_at_read), verbose)
    } else if(ref_class=="ins"){
      cat_at_pos <- 
        read_structure[which(read_structure$map_end==ref_pos)+1,] %>%
        .[which(.$type=="I"),] 
      if(nrow(cat_at_pos)==0){
        mut_at_read <- 0
        base <- NA
        qual <- NA
      } else {  
        length_ins <- str_replace(ref_alt, paste0("^", ref_ref), "") %>%
          nchar()
        if(length_ins==cat_at_pos[1,]$width){
          mut_at_read <- 1
        } else {
          mut_at_read <- -2
        }
        base <- cat_at_pos[1,]$seq
        qual <- cat_at_pos[1,]$qual
      }
    } else if(ref_class=="del"){
      cat_at_pos <- 
        read_structure[which(read_structure$map_end==ref_pos)+1,]%>%
        .[which(.$type=="D"),]
      
      if(nrow(cat_at_pos)==0){
        mut_at_read <- 0
        base <- NA
        qual <- NA
      } else {  
        length_del <- str_replace(ref_ref, paste0("^", ref_alt), "") %>%
          nchar()
        base <- NA
        qual <- NA
        if(length_del==cat_at_pos[1,]$width){
          mut_at_read <- 1
        } else {
          mut_at_read <- -2
        }
      }
    }
  }
  return(tibble(base=base,
                qual=qual,
                mut_at_read=mut_at_read))
}
parse_cigar <- function(bam, qname){
  paired_reads <- bam[which(bam$qname==qname)] %>%
    .[which(as.character(seqnames(.)) %in%
              allowed_inputs("chrom_names"))] %>%
    as_tibble() %>%
    rownames_to_column("mate") %>%
  #parsed_read <- 
    apply(., 1, function(READ){
    make_read_overview_new(as.numeric(READ[["start"]]),
                           READ[["seq"]],
                           READ[["cigar"]],
                           READ[["qual"]]) %>%
      mutate(mate=READ[["mate"]],
             mapq=READ[["mapq"]],
             origin=READ[["origin"]])
  }) %>% bind_rows()  %>% 
    rownames_to_column("id") %>%
    return()
}
evaluate_base <- function(base_info, ref_alt, ref_ref){
  if(base_info$class=="snv"){
    detected <- case_when(
      base_info$base==ref_alt ~ 1,
      base_info$base==ref_ref ~ 0,
      TRUE ~ -2
    )
  } else if(base_info$class=="ins"){
    detected <- case_when(
      base_info$base==ref_alt ~ 1,
      base_info$base==ref_ref ~ 0,
      TRUE ~ -2
    )
  } else {
    detected <- case_when(
      base_info$len_indel==abs(nchar(ref_alt)-nchar(ref_ref)) ~ 1,
      is.na(base_info$len_indel) ~ 0,
      TRUE ~ -2
    )
  }
  return(detected)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom purrr set_names
#' @importFrom dplyr case_when bind_rows as_tibble mutate
#' @importFrom GenomicRanges seqnames
core_tool <- function(qname, bam,
                      ref_pos1, ref_pos2,
                      ref_alt1, ref_alt2,
                      ref_ref1, ref_ref2,
                      ref_class1, ref_class2,
                      verbose=FALSE, version="old"){
  . <- NULL
  ## parse read according to cigar string
  parsed_read <- parse_cigar(bam, qname)
  ## extract base at reference position
  base_info1 <- extract_base_at_refpos(parsed_read, ref_pos1, ref_class1, 
                                           ref_alt1, ref_ref1)
  base_info2 <- extract_base_at_refpos(parsed_read, ref_pos2, ref_class2, 
                                           ref_alt2, ref_ref2)
  ## assign final status to read
  if(is.na(base_info1$base)|is.na(base_info2$base)){
    final_assignment <- "skipped"
  } else {
    mut1_in_read <- evaluate_base(base_info1, ref_alt1, ref_ref1)
    mut2_in_read <- evaluate_base(base_info2, ref_alt2, ref_ref2)
    final_assignment <- case_when(
      sum(mut1_in_read, mut2_in_read)==2 ~ "both",
      sum(mut1_in_read, mut2_in_read)==0 ~ "none",
      #sum(mut1_in_read, mut2_in_read)<(-5) ~ "spanned_out",
      mut1_in_read==1 ~ "mut1",
      mut2_in_read==1 ~ "mut2",
      TRUE ~ "dev_var"
    )    
  }
  return(c(qname=qname, result=final_assignment,
           origin=unique(parsed_read$origin),
           baseq1=base_info1$qual,
           mapq1=base_info1$mapq,
           baseq2=base_info2$qual,
           mapq2=base_info2$mapq
           )
         )
}
#' calculates how many copies are affected by a germnline small variant
#'
#' @param af Allele-frequency of the variant (numeric value between 0 and 1)
#' @param tcn total-copynumber at position of the variant (numeric value >0)
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param chr chromosome of the variant (either format 1,2,..,X,Y or 
#' chr1,..,chrX)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param c_normal expected copy number at position of the variant in normal 
#' tissue, 1 for gonosomes in male samples, and 2 for male autosomes and all 
#' chromosomes in female samples. (The function can also assess the c_normal 
#' parameter by itself, but then the following two inputs must be provided: 
#' chr and sex)
#' @param af_normal Allele-frequency in normal tissue (numeric value between 0 
#' and 1) 0.5 represents heterozygous variants in diploid genome, 
#' 1 would be homozygous. Could be relevant if germline CNVs are present at the 
#' position. Then also the c_normal parameter would have to be adjusted.
#' @return A numeric value indicating the affecting copies for the variant
#' @export
#' @examples 
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#' aff_germ_copies(af=0.67, tcn=2, purity=0.9, chr="chrX", sex="female")
aff_germ_copies <- function(chr, af, tcn, purity, sex, 
                            c_normal=NULL, af_normal=0.5){
  cv <- formula_checks(chr, af, tcn, purity, sex, c_normal, af_normal)
  aff_cp <- cv$af*cv$tcn+(cv$af-cv$af_normal)*cv$c_normal*((1/cv$purity)-1)
  return(aff_cp) 
  # alternative: (af*(purity*tcn+c_normal*(1-purity))-cc*(1-purity))/purity 
}
#' calculates how many copies are affected by a somatic small variant
#'
#' @param af Allele-frequency of the variant (numeric value between 0 and 1)
#' @param tcn total-copynumber at position of the variant (numeric value >0)
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param chr chromosome of the variant (either format 1,2,..,X,Y or 
#' chr1,..,chrX)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param c_normal expected copy number at the position of the variant in 
#' normal tissue, 1 for gonosomes in male samples, and 2 for male autosomes and 
#' all chromosomes in female samples. (The function can also assess the 
#' c_normal parameter by itself, but then the following two inputs must be 
#' provided: chr and sex)
#' @return A numeric value indicating the affecting copies for the variant
#' @examples 
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#' aff_som_copies(chr="chrX", af=0.67, tcn=2, purity=0.9, sex="female")
#' @export
aff_som_copies <- function(chr, af, tcn, purity, sex, c_normal=NULL){
  cv <- formula_checks(chr, af, tcn, purity, sex, c_normal)
  aff_cp <- cv$af*(cv$tcn+cv$c_normal*((1/cv$purity)-1))
  return(aff_cp)
  ## alternative: aff_cp = af*(tcn+c_normal/purity-c_normal)
}
#' @keywords internal
#' description follows
formula_checks <- function(chr, af, tcn, purity, sex, c_normal, af_normal=0.5){
  purity <- check_purity(purity)
  af <- check_af(af)
  af_normal <- check_af(af_normal)
  tcn <- check_tcn(tcn)
  if(is.null(c_normal)){
    sex <- check_sex(sex)
    chr <- check_chr(chr)
    if((sex=="male"&(chr=="X"|chr=="chrX"))|(chr=="Y"|chr=="chrY")){
      c_normal <- 1
    } else {
      c_normal <- 2
    }
  } else {
    c_normal <- check_ploidy(c_normal)
  }
  return(list(af=af, tcn=tcn, purity=purity, c_normal=c_normal, 
              af_normal=af_normal))
}
#' @keywords internal
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate select mutate_at rowwise
#' @importFrom stringr str_detect
merge_sCNAs <- function(obj, somCna){
  alt <- af <- tcn <- tcn_assumed <- cna_type <- gene <- ref <- NULL
  obj %>%
    mergeByOverlaps(somCna) %>% as_tibble() %>%
    mutate(cna_type=ifelse(str_detect(cna_type, "LOH"), "LOH", "HZ")) %>%
    select(chr=1, pos=2, gene, ref, alt, af, tcn, cna_type, tcn_assumed) %>%
    mutate_at(.vars = c("af", "tcn"), .funs=as.numeric) %>%
    rowwise() %>%
    return()
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate filter
prepare_somatic_variant_table <- function(somSmallVars, templateGenes, 
                                          somCna, purity, sex){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- 
    aff_cp <- wt_cp <- . <- tcn_assumed <- NULL
  if(!is.null(somSmallVars)){
    rel <- somSmallVars%>%
      merge_sCNAs(., somCna) %>%
      
      ################################
      mutate(
        origin="somatic",
        class=define_class(ref, alt),
        aff_cp = aff_som_copies(chr, af, tcn, purity, sex),
        wt_cp = tcn-aff_cp,
        pre_info = 
          ifelse(isTRUE(str_detect(cna_type,'LOH')),
            ifelse(wt_cp<=0.5,
              paste(
                  'somatic-variant -> LOH -> left wt copies',round(wt_cp,2),
                  '-> all copies affected'),
              paste(
                  'somatic-variant -> LOH -> ',
                  "left wt copies",round(wt_cp,2),"-> wt copies left")),
            ifelse(chr=="X"&sex=="male",
              ifelse(wt_cp<=0.5,
                paste(
                  "somatic-variant -> chrX & male -> ",
                  "left wt copies",round(wt_cp,2),"-> all copies affected"),
                paste(
                  'somatic-variant -> chrX & male -> left wt copies ',
                  round(wt_cp,2),
                  '-> wt copies left')),
                'somatic-variant -> no LOH -> wt copies left')
        )
      )%>%
      filter(gene %in% templateGenes)
    if(nrow(rel)==0){
      return(NULL)
    } else {
      return(rel)
    }
  } else {
    return(NULL)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr compact
#' @importFrom dplyr filter as_tibble select bind_rows mutate
extract_all_dels_of_sample <- function(somCna, geneModel, DEL_TYPE, 
                                       byTcn, sex, ploidy, include_dels){
  tcn <- cna_type <- seqnames <- gene <- . <- NULL
  if(include_dels==TRUE){
    TBL_CNV <- somCna %>%
      as_tibble() %>%
      mutate(tcn=as.numeric(tcn))
    if(DEL_TYPE=="homdel"){
      if(byTcn==TRUE){
        relevant_CNAs <- TBL_CNV %>%
          filter(round(tcn)==0)
      } else {
        relevant_CNAs <- TBL_CNV %>%
          filter(str_detect(cna_type, allowed_inputs("cna_homdel_annotation"))) 
      }
      pre_info <- 'homdel -> all copies affected'
    } else {
      if(byTcn==TRUE){
        if(sex=="female"){
          relevant_CNAs <- TBL_CNV %>%
            filter(round(tcn)<ploidy)
        } else {
          relevant_CNAs <- TBL_CNV %>%
            filter(!str_detect(seqnames, "X|Y")&round(tcn)<ploidy)
        }
      } else {
        relevant_CNAs <- TBL_CNV %>%
          filter(!str_detect(cna_type, 
                             allowed_inputs("cna_homdel_annotation"))) %>%
          filter(str_detect(cna_type, 
                            allowed_inputs("cna_incompletedel_annotation")))
      }
      pre_info <- "incompletedel -> wt-copies left"
    } 
    if(nrow(relevant_CNAs)==0){
      return(NULL)
    } else {
      merged_CNAs <- apply(relevant_CNAs, 1, function(CNV){
        GR_single_CNV <- GRanges(as_tibble(t(CNV)))
        suppressWarnings(
        rel <- subsetByOverlaps(geneModel, GR_single_CNV) %>%
          as_tibble()
        )
        if(nrow(rel)==0){
          return(NULL)
        } else {
          rel %>%
            select(gene) %>%
            mutate(chr=CNV[["seqnames"]],
                   pos=as.numeric(CNV[["start"]]),
                   tcn=as.numeric(CNV[["tcn"]]),
                   tcn_assumed=FALSE,
                   af=NA,
                   origin="somatic",
                   class=DEL_TYPE,
                   ref= NA,
                   alt= NA,
                   cna_type=DEL_TYPE,
                   aff_cp=ploidy-as.numeric(CNV[["tcn"]]),
                   wt_cp=as.numeric(CNV[["tcn"]]),
                   pre_info=pre_info,
            ) %>%
            return()
        }
      })%>% compact() %>% bind_rows() %>% unique()
      if(nrow(merged_CNAs)==0){
        return(NULL)
      } else {
        return(merged_CNAs)
      }
    }
  } else {
    return(NULL)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom dplyr mutate select
prepare_germline_variants <- function(germSmallVars, somCna, purity, sex){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- aff_cp <- 
    origin <- pos <- wt_cp <- pre_info <- . <- tcn_assumed <- NULL
  if(is.null(germSmallVars)){
    return(NULL)
  } else {
    df_germ <- germSmallVars %>% 
      merge_sCNAs(., somCna) %>%
      mutate(
        origin="germline",
        class=define_class(ref, alt),
        aff_cp = aff_germ_copies(chr, af, tcn, purity, sex, 1),
        wt_cp = tcn-aff_cp,
        pre_info = 
          ifelse(isTRUE(cna_type=='LOH'),
            ifelse(af>=0.5,
                 paste(
                   'germline-variant -> LOH ->',
                   'AF >= 0.5 -> all copies affected'),
                 paste(
                   'germline-variant -> LOH ->',
                   ' AF < 0.5 -> variant lost in tumor')),
            'germline-variant -> no LOH -> wt copies left')
      ) %>%
      select(gene,
             origin,
             class,
             chr,
             pos,
             ref,
             alt,
             af,
             tcn, 
             tcn_assumed,
             cna_type,
             aff_cp,
             wt_cp,
             pre_info
      ) 
    return(df_germ)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_match
#' @importFrom GenomicRanges elementMetadata elementMetadata<-
assign_correct_colnames <- function(obj, type){
  . <- NULL
  if(type=="scna"){
    col_tcn <- str_match(
      nm_md(obj), paste(allowed_inputs("colnames_tcn"), collapse="|")) %>%
      .[which(!is.na(.))]
    elementMetadata(obj)[,"tcn"] <- 
      elementMetadata(obj)[,col_tcn]   
    if("LOH" %in% nm_md(obj)){
      obj$cna_type <- case_when(obj$LOH==TRUE ~ "LOH",
                                    TRUE ~ "HZ")
    } else {
      col_cna_type <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_cna_type"), 
                          collapse="|")) %>%
        .[which(!is.na(.))]   
      elementMetadata(obj)[,"cna_type"] <- 
        elementMetadata(obj)[,col_cna_type]
    }
  } else {
    col_gene <- str_match(
      nm_md(obj), 
      paste(allowed_inputs("colnames_gene"), collapse="|")) %>%
      .[which(!is.na(.))]
    elementMetadata(obj)[,"gene"] <- 
      elementMetadata(obj)[,col_gene]  
    if(type=="small_vars"){
      col_af <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_af"), collapse="|")) %>%
        .[which(!is.na(.))]
      col_ref <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_ref"), collapse="|")) %>%
        .[which(!is.na(.))]
      col_alt <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_alt"), collapse="|")) %>%
        .[which(!is.na(.))]
      elementMetadata(obj)[,"af"] <- 
        elementMetadata(obj)[,col_af]
      elementMetadata(obj)[,"ref"] <- 
        elementMetadata(obj)[,col_ref]
      elementMetadata(obj)[,"alt"] <- 
        elementMetadata(obj)[,col_alt]
      elementMetadata(obj)["mid"] <- 
        seq_len(length(obj))    
    }    
  }

  return(obj)
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_split
#' @importFrom purrr map_chr set_names
make_phasing_combinations <- function(df_mut_to_combine){
  . <- pos1 <- pos2 <- mut_id1 <- mut_id2 <- NULL
  outer(df_mut_to_combine$pos, df_mut_to_combine$pos, `paste`) %>%
    .[which(upper.tri(.))] %>%
    as.data.frame() %>% 
    set_names(nm='raw') %>%
    mutate(
      pos1=str_split(raw, ' ') %>%map_chr(.,1) %>% as.numeric(),
      pos2=str_split(raw, ' ') %>%map_chr(.,2) %>% as.numeric(),
      dist = abs(pos1-pos2)
    ) %>%
    # add all information about every mutation to the data frame
    merge(df_mut_to_combine %>% 
            set_names(nm=names(df_mut_to_combine) %>% paste0(.,'1')), 
          by='pos1', all.x=TRUE) %>%
    merge(df_mut_to_combine %>% 
            set_names(nm=names(df_mut_to_combine) %>% paste0(.,'2')), 
          by='pos2', all.x=TRUE) %>%
    .[order(.$dist),] %>%
    rowwise() %>%
    mutate(
      comb_id=paste(mut_id1, mut_id2, sep='-'),
      comb_id_sorted=paste(paste0("m",
                           sort(as.numeric(str_match(
                             c(mut_id1, mut_id2), "\\d+")))),
                    collapse = "-")
      ) %>%
    ungroup()
}
#' @keywords internal
#' description follows
check_for_overlapping_reads <- function(bamDna, bamRna,
                                        ref_chr1, 
                                        ref_chr2, 
                                        ref_pos1, 
                                        ref_pos2, 
                                        verbose=FALSE){
  #vm("checking DNA", verbose)
  dna_bam <- prepare_raw_bam_file(bamDna, 
                                  ref_chr1, 
                                  ref_chr2, 
                                  ref_pos1, 
                                  ref_pos2,
                                  verbose) 
  if(length(dna_bam)!=0){
    dna_bam$origin <- "DNA"
  }
  if(!is.null(bamRna)){
    #vm("checking RNA", verbose)
    rna_bam <- prepare_raw_bam_file(bamRna, 
                                    ref_chr1, 
                                    ref_chr2, 
                                    ref_pos1, 
                                    ref_pos2,
                                    verbose)
    if(length(rna_bam)==0){
      rna_bam <- NULL
    } else {
      rna_bam$origin <- "RNA"
    }
  } else {
    rna_bam <- NULL
  }
  return(c(dna_bam, rna_bam))
}
ascii_to_dec <- function(ascii_encoded){
  ## apply statistics
  #ascii_encoded <- str_sub(qual_string, start=1, end=1)
  exp <- (as.integer(charToRaw(ascii_encoded))-33)/(-10)
  # Convert the QUAL string to raw hexadecimal
  return(10^exp)
}

#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr tibble 
classify_reads <- function(line, bamDna, bamRna, verbose){
  dist <- abs(as.numeric(line[['pos1']])-as.numeric(line[['pos2']]))
  ref_pos1 <- as.numeric(line[['pos1']])
  ref_pos2 <- as.numeric(line[['pos2']])
  ref_chr1 <- as.character(line[['chr1']])
  ref_chr2 <- as.character(line[['chr2']])
  ref_alt1 <- as.character(line[['alt1']])
  ref_alt2 <- as.character(line[['alt2']])
  ref_ref1 <- as.character(line[['ref1']])
  ref_ref2 <- as.character(line[['ref2']])
  ref_class1 <- as.character(line[['class1']])
  ref_class2 <- as.character(line[['class2']])
  #vm(line, verbose)
  bam <- check_for_overlapping_reads(bamDna,
                                     bamRna,
                                     ref_chr1,
                                     ref_chr2,
                                     ref_pos1,
                                     ref_pos2, 
                                     verbose)
  #print(bam)
  #print(ref_ref1)
  #print(ref_ref2)
  if(!length(bam)==0){  
    vm("reads detected, applying core function", verbose)
    classified_reads <- lapply(unique(bam$qname),
                               core_tool,
                               bam, 
                               ref_pos1,
                               ref_pos2,
                               ref_alt1,
                               ref_alt2,
                               ref_ref1,
                               ref_ref2,
                               ref_class1,
                               ref_class2,
                               verbose) %>%
      bind_rows() #%>%
    #  rowwise() %>%
    #  mutate(baseq1_conv=ascii_to_dec(baseq1),
    #         baseq2_conv=ascii_to_dec(baseq2))
    
    
    
  } else {
    #vm("no reads detected, return empty df", verbose)
    classified_reads <- tibble()
  } 
  return(classified_reads)
}
load_vars_from_vcf <- function(vcf, df_gene, gr_roi){
  vars_in_between_raw <- lapply(vcf, function(VCF){
    tab_vcf <- TabixFile(VCF)
    vcf_region <- 
      readVcf(tab_vcf, "hg19", 
              param=gr_roi) %>%
      rowRanges() %>%
      return()
  }) %>%
    Reduce(function(x,y)c(x,y),.) %>%
    as_tibble() %>%
    ## maybe improve handling of more entries in DNAstringLists
    ## now it takes the first
    mutate(ALT=unlist(lapply(ALT, function(x){as.character(x[[1]][1])}))) %>%
    rowwise() %>%
    mutate(ALT=as.character(ALT),
           class=define_class(REF, ALT)) %>%
    select(1,2,3, ref=REF, alt=ALT, class) %>%
    unique() %>%
    select(-end) %>%
    dplyr::rename(pos=start, chr=seqnames) %>%
    left_join(df_gene %>% select(chr, pos, mut_id),
              by=c("chr"="chr",
                   "pos"="pos")) %>%
    arrange(mut_id) %>%
    ungroup() 
  return(vars_in_between_raw)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom Rsamtools TabixFile
#' @importFrom VariantAnnotation readVcf
#' @importFrom DelayedArray rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr rename arrange as_tibble left_join mutate ungroup rowwise ungroup select
check_for_snps_between_main_muts <- function(main_comb, vcf, df_gene,
                                             snp_dist=0){
  . <- ALT <- REF <- end <- start <- seqnames <- chr <- pos <- mut_id <- NULL
  pos_minor <- as.numeric(min(main_comb[['pos2']],
                   main_comb[['pos1']]))-snp_dist
  pos_major <- as.numeric(max(main_comb[['pos2']],
                   main_comb[['pos1']]))+snp_dist
  gr_roi <- GRanges(paste0(main_comb[["chr1"]], ":",
                                          pos_minor, "-", pos_major))
  vars_in_between_raw <- load_vars_from_vcf(vcf, df_gene, gr_roi) %>%
 
    ## remove main mutations if they are in between... because would be
    ## duplicated with upper analysis
    filter(is.na(mut_id)|
             mut_id %in% c(main_comb[["mut_id1"]], main_comb[["mut_id2"]])) %>%
    mutate(
      pos = as.numeric(pos),
      site=case_when(
        between(pos, pos_minor, pos_major) ~ "between",
          pos>pos_major ~ "downstream",
          pos<pos_minor ~ "upstream"
    ))
  return(vars_in_between_raw)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr mutate as_tibble
check_read_presence <- function(sub_comb, bamDna, bamRna){
  bam_raw <- check_for_overlapping_reads(bamDna, bamRna,
                                         sub_comb[["chr1"]],
                                         sub_comb[["chr2"]],
                                         sub_comb[["pos1"]],
                                         sub_comb[["pos2"]]) 
  if(length(bam_raw)!=0){
    bam <- bam_raw 
    bam$comb_id <- sub_comb[["comb_id"]]
  } else {
    bam <- NULL
  }
  res_tibble <- as_tibble(t(sub_comb)) %>%
    mutate(nreads=0.5*length(bam_raw))
  rl <- list(rt=res_tibble, bam=bam)
  return(rl)
}
#' @keywords internal
#' @importFrom igraph graph_from_data_frame shortest_paths
find_path <- function(connections, mut_id1, mut_id2) {
  graph <- graph_from_data_frame(connections, directed = FALSE)
  shortest_path <- 
    shortest_paths(graph, from = mut_id1, to = mut_id2, mode = "all")$vpath
  if (length(shortest_path) > 0) {
    return(names(unlist(shortest_path)))
  } else {
    return(NULL)
  }
}

#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr mutate
return_null_result <- function(main_comb, bamRna, ext_snp_info){
  RESULT <- rep(0, 9) %>% t() %>% as.data.frame() %>% 
    set_names(c("both", "none", "mut1", "mut2", "dev_var", 
                "no_overlap", "none_raw", "DNA_rds", "RNA_rds")) %>% 
    mutate(
      comb=main_comb[['comb_id']],
      class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                       sep="-"),
      dist=as.numeric(main_comb[["dist"]]),
      status="null",
      score=1,
      wt_cp=min(as.numeric(main_comb[['tcn1']]), 
                as.numeric(main_comb[['tcn2']])) %>% 
        round(2),
      tcn=min(as.numeric(main_comb[['tcn1']]), 
              as.numeric(main_comb[['tcn2']])) %>% 
        round(2),
      info=ifelse(!is.null(bamRna),
                  paste('unclear: no ovrlp. reads;', ext_snp_info),
                  paste('unclear: no ovrlp. reads & no RNA file;', 
                        ext_snp_info))
    )
} 

get_genotype <- function(gt, status){
  if(status=="same"){
    return(gt)
  } else if(gt=="1|0"){
    return("0|1")
  } else if(gt=="0|1"){
    return("1|0")
  } 
}


perform_level_2_phasing <- function(df_gene, vcf, haploBlocks, phasedVcf,
                         bamDna, bamRna, 
                        purity, snp_dist=10000, distCutOff,
                        printLog, verbose,
                        geneDir){
  
  catt(printLog, 4, "Initializing level 2 phasing")
  #print(paste("GENE:", unique(df_gene$gene)))
  ## annotate to df_gene in which haploblock the variant is located, 0 means no 
  ## haploblock

  
  df_gene_hap <- apply(df_gene, 1, function(mut){
    gr_roi <- GenomicRanges::GRanges(paste0(mut[["chr"]], ":", 
                                              as.numeric(mut[["pos"]]), "-",
                                                as.numeric(mut[["pos"]])))
    #print(gr_roi)
    #print(haploBlocks)
    hap <- IRanges::subsetByOverlaps(haploBlocks, gr_roi)
    if(length(hap)==0){
      det_hap_id <- 0
    } else {
      det_hap_id <- hap$hap_id
      #print(mut)
      #print(hap)
    }
    return(as_tibble(t(mut)) %>% mutate(hap_id=det_hap_id))
  }) %>% bind_rows()
  
  #print(df_gene_hap)
  ## check if at least two variants are in the same haploblock
  if(any(df_gene_hap$hap_id > 0)&any(table(
    df_gene_hap[which(df_gene_hap$hap_id>0),]$hap_id
    )>1)){
    vm("at least two variants detected in same haploblock", verbose)
    #write_tsv(df_gene_hap, file=file.path("/home/m168r/home_extension/ZP_revision", paste0(unique(df_gene$gene), ".tsv")))
    relevant_haploblocks <- table(df_gene_hap$hap_id) %>% .[which(.>1)] %>% names() %>% .[which(!.=="0")]
    #print("relevant haploblocks")
    #print(relevant_haploblocks)
    ## iterate over haploblocks
    per_haploblock <- lapply(relevant_haploblocks, function(HAP_ID){
      #print(haploBlocks)
      #print(HAP_ID)
      region_to_load <- haploBlocks[which(haploBlocks$hap_id==HAP_ID)]
      
      vcf_to_load <- phasedVcf %>% .[which(str_detect(., paste0("chr",unique(df_gene$chr), ".vcf")))]
      #print(region_to_load)
      #print("loading vcf")
      #print(vcf_to_load)
      vcf_region <- 
        VariantAnnotation::readVcf(vcf_to_load, "hg19")
      
      #print("vcf loaded")
      gt <- VariantAnnotation::geno(vcf_region)$GT %>% as.character()
      comb_vcf <- rowRanges(vcf_region)
      comb_vcf$gt <- gt
      
      #print(comb_vcf)
      #print(region_to_load)
      
      phased_snps_in_haploblock <- subsetByOverlaps(comb_vcf, region_to_load) %>%
        as_tibble() %>%
        filter(str_detect(gt, "\\|")) %>%
        filter(!gt=="1|1") 
      
      if(nrow(phased_snps_in_haploblock)==0){
        vm("no phased snps/heterozygous SNPS in haploblock", verbose)
        return(NULL)
      } else {
        phased_snps_in_haploblock_ids <- phased_snps_in_haploblock %>%
          mutate(snp_id=paste0("s", c(1:nrow(.))))
        if(!is.null(geneDir)){
          write_tsv(phased_snps_in_haploblock_ids, 
                    file=file.path(geneDir, 
                                   paste0("all_snps_hb_", 
                                          HAP_ID, ".tsv")))
        }
        #print(phased_snps_in_haploblock)
        
        ## get every mutation located in the haploblock
        muts_in_hap <- df_gene_hap %>%
          filter(hap_id==HAP_ID)
        
        ## calculate distance of each snp to each variant
        df_dist <- apply(muts_in_hap, 1, function(mut){
          
          phased_snps_in_haploblock_ids %>%
            mutate(dist=abs(start-as.numeric(mut[["pos"]])),
                   mut_id=mut[["mut_id"]]) %>%
            return()
          
        }) %>% bind_rows() %>%
          filter(dist<distCutOff) %>%
          arrange(dist) %>%
          mutate(comb_id=paste(mut_id, snp_id, sep="-"))
        
        #print(df_dist)
        if(!is.null(geneDir)){
          write_tsv(df_dist, 
                  file=file.path(geneDir, 
                                 paste0("full_dist_df_hb_", 
                                        HAP_ID, ".tsv")))
        }
        if(length(unique(df_dist$mut_id))>1){
          phased_per_mut_full <- lapply(unique(df_dist$mut_id), function(MUT_ID){
            #message("phasing", MUT_ID)
            df_dist_mut <- df_dist %>%
              filter(mut_id==MUT_ID)
            #print(df_dist_mut)
            high_conf_result <- FALSE
            i <- 1
            phased_muts <- tibble()
            ## only for documentation
            sub_status_list <- tibble()
            while(high_conf_result==FALSE&i<=nrow(df_dist_mut)){
              cat(i)
              
              comb <- df_dist_mut[i,]
              #print("hier kommt comb")
              #print(comb)
              sub_comb <- bind_rows(
                df_gene %>% 
                  filter(mut_id==comb$mut_id) %>%
                  select(chr, pos, mut_id, ref, alt, class),
                comb %>%
                  mutate(ALT=unlist(
                    lapply(ALT, function(x){as.character(x[[1]][1])})),
                    class="snv") %>%
                  select(chr=seqnames, pos=start, mut_id=snp_id, ref=REF, alt=ALT, 
                         class)
              ) %>%
                make_phasing_combinations()
              
              ## the function classify_reads is usually used inside of apply()
              ## thats why the input needs to be adjusted in the format of apply iterations:
              comb_vec <- as.character(sub_comb) %>% set_names(nm=names(sub_comb))
              
       
              
              res <- classify_reads(comb_vec, bamDna, bamRna, verbose)
              
              if(nrow(res)!=0){
                sub_status_combination <- 
                  classify_combination(res,
                                       purity,
                                       eval_full=FALSE,
                                       printLog)
                #vm("111111111111", verbose)
                sub_status_list <- bind_rows(sub_status_list,
                                             sub_status_combination %>% 
                                               mutate(comb_id=comb$comb_id))
                #print(sub_status_list)
                #vm("222222221", verbose)
                if(!sub_status_combination$status=="null"){
                  phased <- tibble(mut_id=comb$mut_id,
                                   gt=get_genotype(comb$gt, sub_status_combination$status))
                  phased_muts <- bind_rows(phased_muts, phased)
                  
                  high_conf_result <- TRUE
                  
                } 
                
              }
              i <- i +1 

            } 

            
            return(list(phased_muts, sub_status_list))
            
            
            # if(length(bam_raw)>0){
            #   write_tsv(tibble(x="hello"), file=file.path("/home/m168r/home_extension/ZP_revision", paste0(unique(df_gene$gene), ".tsv")))
            # }
            
          }) #%>% bind_rows()
          #vm("333333333", verbose)
          phased_per_mut <- bind_rows(lapply(phased_per_mut_full, nth, 1))
          
          #print("done phasing of variants to snps")
          #print(phased_per_mut)
          if(nrow(phased_per_mut)>0){
             #print(phased_per_mut_full)
            if(!is.null(geneDir)){
                log_phasing_results <- bind_rows(lapply(phased_per_mut_full, nth, 2))
                write_tsv(log_phasing_results, 
                          file=file.path(geneDir, 
                                         paste0("phasing_info_hb_", 
                                                HAP_ID, ".tsv")))
            }
            
            #vm("111111111111", verbose)
            status_per_comb <- outer(phased_per_mut$mut_id, 
                                     phased_per_mut$mut_id, `paste`) %>%
              .[which(upper.tri(.))] %>%
              as.data.frame() %>% 
              set_names(nm='raw') %>%
              mutate(
                mut_id1=str_split(raw, " ") %>% map_chr(.,1),
                mut_id2=str_split(raw, " ") %>% map_chr(.,2)
              ) %>%
              left_join(phased_per_mut %>% select(mut_id1=mut_id, gt1=gt), 
                        by="mut_id1") %>%
              left_join(phased_per_mut %>% select(mut_id2=mut_id, gt2=gt), 
                        by="mut_id2") %>%
              rowwise() %>%
              mutate(status=case_when(gt1==gt2 ~ "same",
                                      TRUE ~ "diff"),
                     #comb_id=paste(mut_id1, mut_id2, sep='-')
                     comb_id=paste(paste0("m",
                                          sort(as.numeric(str_match(
                                            c(mut_id1, mut_id2), "\\d+")))),
                                   collapse = "-")
                     ) %>%
              ungroup() %>%
              select(comb_id, status)  
          } else {
            status_per_comb <- NULL
          }
         
          #vm("222222222", verbose)
          return(status_per_comb)
          #print(table(phased_per_mut))
          
          
          # if(nrow(phased_per_mut)>1){
          #   write_tsv(tibble(x="hello"), file=file.path("/home/m168r/home_extension/ZP_revision/gene_out", paste0(unique(df_gene$gene), ".tsv")))
          # }        
        } else {
          vm("only one mut closer than dist limit to snp", verbose)
          return(NULL)
        }
      }
    })  %>%
      bind_rows()

    if(nrow(per_haploblock)==0){
      return(NULL)
    } else {
      return(per_haploblock)
    }
  } else {
    #print("no pair of variants located in the same haploblock")
    return(NULL)
  }
}

perform_direct_phasing <- function(all_combinations, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail){
  final_combinations <-  apply(all_combinations, 1, function(main_comb){
    
    catt(printLog, 2 ,c("combination:", main_comb[["comb_id"]]))
    return_null_result <- TRUE
    #try_snp_phasing <- TRUE
    level_2_phasing <- TRUE
    #man_snp_phasing <- FALSE
    #ext_snp_phasing <- NULL
    RESULT <- NULL
    full_read_info <- NULL
    #vm("classifying reads for main mutations", verbose)
    main_classified_reads <- classify_reads(main_comb, bamDna, bamRna, verbose)
    catt(printLog, 3, 
         c(nrow(main_classified_reads), 
           "reads / read-pairs overlapping both positions"))
    if(nrow(main_classified_reads)!=0){
      return_null_result <- FALSE
      if(showReadDetail==TRUE){
        full_read_info <- main_classified_reads %>%
          select(qname, result, origin) %>%
          mutate(comb=main_comb[['comb_id']])
      } else {
        full_read_info <- NULL
      }
      RESULT <- classify_combination(main_classified_reads,
                                     purity,
                                     eval_full = TRUE, 
                                     printLog,
                                     main_comb[['tcn1']],
                                     main_comb[['tcn2']],
                                     main_comb[['aff_cp1']],
                                     main_comb[['aff_cp2']]
      ) %>%  
        mutate(
          DNA_rds=nrow(main_classified_reads %>% filter(origin=="DNA")),
          RNA_rds=nrow(main_classified_reads %>% filter(origin=="RNA")),
          dist=as.numeric(main_comb[["dist"]]),
          class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                           sep="-"),
          comb=main_comb[['comb_id']]
        ) 
      if(RESULT$status!="null"){
        #try_snp_phasing <- FALSE
        level_2_phasing <- FALSE
        return_null_result <- FALSE
        catt(printLog, 3, c("status:", RESULT$status))
      } else {
        catt(printLog, 3, "status can not be defined")
      }
    }   
    # if(is.null(vcf)){
    #   ext_snp_info <- "no vcf provided for extended SNP phasing"
    # } else if(try_snp_phasing==FALSE){
    #   ext_snp_info <- ""
    # } else if(man_snp_phasing==TRUE&as.numeric(main_comb[["dist"]])<distCutOff){
    #   ext_snp_raw <- perform_extended_snp_phasing(main_comb, vcf,
    #                                               df_gene, bamDna, bamRna, 
    #                                               purity, full_read_info, RESULT,
    #                                               ext_snp_phasing,
    #                                               printLog, verbose)
    #   ext_snp_phasing <- ext_snp_raw[["ext_snp_phasing"]]
    #   ext_snp_info <- ext_snp_raw[["ext_snp_info"]]
    #   RESULT <-  ext_snp_raw[["RESULT"]]
    # } else { 
    #   ext_snp_info <- "muts too far away for SNP phasing"
    # }
    # catt(printLog, 3, ext_snp_info)
    if(return_null_result==TRUE){
      RESULT <- return_null_result(main_comb, bamRna, "test")
      full_read_info <- NULL
    }
    return(list(assigned=RESULT, detailed=full_read_info, 
                #snp_phasing=ext_snp_phasing, 
                level_2_phasing=level_2_phasing))
  })
  return(final_combinations)
}
#' @keywords internal
#' function to decide if indirect phasing should and can be performed
eval_direct_results <- function(status_per_comb, phasedVcf,
                    haploBlocks, mode="fast"){
  try <- FALSE
  if(!is.null(phasedVcf)&!is.null(haploBlocks)){
    if(mode=="fast"){
      if("null" %in% status_per_comb$status){
        if(!("diff" %in% status_per_comb$status)|
           ("same" %in% status_per_comb$status)){
          try <- TRUE
        }
      }      
    } 
  }
  return(try)
}

select_status <- function(status_direct, status_indirect){
  if("diff" %in% c(status_direct, status_indirect)){
    return("diff")
  } else if("same" %in% c(status_direct, status_indirect)){
    return("same")
  } else {
    return("null")
  }
}

#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
phase <- function(df_gene, bamDna, bamRna, 
                  showReadDetail, purity, vcf, haploBlocks, phasedVcf,
                  distCutOff, printLog, verbose, logDir){
  pos1 <- pos2 <- mut_id1 <- mut_id2 <- qname.first <- result <- . <- origin <- 
    qname <- mut_id <- nreads <- status <- geneDir <- NULL
  if(!is.null(logDir)){
    geneDir <- file.path(logDir, unique(df_gene$gene))
    dir.create(geneDir)
  }
  #print("DF_GENE:")
  #print(df_gene)
  ## (1): define all combinations of variants to be phased
  all_combinations <- make_phasing_combinations(df_gene) 
  ## (2): perform direct phasing between variants (read-based)
  direct_mut_phasing <- perform_direct_phasing(all_combinations, bamDna, bamRna, 
                                               purity,
                                               verbose, printLog, 
                                               showReadDetail)
  
  status_per_comb <- lapply(direct_mut_phasing, nth, 1) %>% bind_rows() %>%
    select(comb_id=comb, status)
  
  ## (3): check if indirect phasing can/should be performed
  do_indirect_phasing <- eval_direct_results(status_per_comb, phasedVcf,
                                             haploBlocks)
  ## just to pre define... will be overwritten if necesarry
  merged_comb_indirect <- all_combinations %>%
    mutate(status_indirect="null")
  #print(do_indirect_phasing)
  
  if(do_indirect_phasing){
    #print(df_gene)
    indirect_phasing_result <- perform_level_2_phasing(df_gene, vcf, 
                                                      haploBlocks, phasedVcf,
                                                      bamDna, bamRna,
                                                      purity, snp_dist=10000, 
                                                      distCutOff,
                                                      printLog, verbose,
                                                      geneDir)
    vm("indirect phasing done", verbose)
    if(!is.null(indirect_phasing_result)){
      merged_comb_indirect <- left_join(all_combinations,
                                        indirect_phasing_result %>% 
                                          select(comb_id, status_indirect=status),
                                        by=c("comb_id_sorted"="comb_id"))      
    }
  }
  merged_comb_full <- left_join(merged_comb_indirect,
                                status_per_comb %>% 
                                  select(comb_id, status_direct=status),
                           by=c("comb_id"="comb_id")) %>%
    rowwise() %>%
    mutate(status=case_when(status_direct!="null" ~ status_direct,
                            status_indirect!="null" ~ status_indirect,
                            TRUE ~ "null"
                            ),
             #select_status(status_direct, status_indirect),
           tcn=min(tcn1, tcn2),
           wt_cp=calc_left_wt_copies(tcn,
                                     status,
                                              aff_cp1,
                                              aff_cp2),
           score = ifelse(status=="diff"&wt_cp<0.5,
                               2, 1),
           comb=comb_id,
           phased=case_when(status_direct!="null" ~ "direct",
                            status_indirect!="null" ~ "indirect",
                            TRUE ~ "not"
           )
           )
  
  #print(direct_mut_phasing)
  full_phasing_result <- list(
    merged_comb_full,
    direct_mut_phasing
  )
  if(!is.null(geneDir)){
    write_tsv(lapply(direct_mut_phasing, nth, 1) %>% bind_rows(), 
              file=file.path(geneDir, "direct_phasing.tsv"))
  }
  return(full_phasing_result) 
}
eval_phasing <- function(all_comb, df_gene, log_obj, printLog, verbose){
  vm("start evaluation of phasing results", verbose)
  if(nrow(all_comb)==1){
    vm("only one combination", verbose)
    log_obj[['score']] <- all_comb$score
    log_obj[['info']] <- all_comb$info
    log_obj[['tcn']] <- all_comb$tcn
    log_obj[['aff_cp']] <- as.numeric(all_comb$tcn-all_comb$wt_cp)
    log_obj[['class']] <- paste("comb:", all_comb$comb)
  } else if(max(all_comb$score)==2){
    vm("at least one combination that affects all copies", verbose)
    most_important_comb <- all_comb %>% 
      filter(wt_cp==min(all_comb$wt_cp)) %>%
      .[1,]
    log_obj[['score']] <- 2
    log_obj[['info']] <- 
      "at least one of the combinations of muts affects all copies"
    log_obj[['tcn']] <- most_important_comb$tcn
    log_obj[['aff_cp']] <- 
      most_important_comb$tcn-most_important_comb$wt_cp %>% 
      as.numeric()
    log_obj[['class']] <- paste("comb:", most_important_comb$comb)
  } else if(str_count(paste(all_comb$status, collapse=' '),
                      'diff')==nrow(all_comb)&
            nrow(all_comb)>=mean(as.numeric(df_gene$tcn))){
    vm("unclear result: tcn 2, all diff but left_wt copies too high", verbose)
    most_important_comb <- all_comb %>% 
      filter(wt_cp==min(all_comb$wt_cp)) %>%
      .[1,]
    log_obj[['score']] <- 1
    log_obj[['info']] <- paste(
      'unclear: all affected is very likely because all',nrow(all_comb),
      'mutations are at different copies and tcn=', 
      mean(as.numeric(df_gene$tcn)),
      'tumor might be subclonal and therefore wt-copies left')
    log_obj[['tcn']] <- most_important_comb$tcn
    log_obj[['aff_cp']] <- 
      most_important_comb$tcn-most_important_comb$wt_cp %>% 
      as.numeric()
    log_obj[['class']] <- paste("comb:", most_important_comb$comb)
#  } else if(str_detect(paste(all_comb$status, collapse=' '),
    #                    'unclear: muts at different copies and tcn=2')){
    # most_important_comb <- all_comb %>% 
    #   filter(str_detect(status, 
    #                     'unclear: muts at different copies and tcn=2')) %>%
    #   filter(wt_cp==min(all_comb$wt_cp)) %>%
    #   filter(rownames(.)==1)
    # log_obj[['score']] <- 1
    # log_obj[['info']] <- paste(
    #   "unclear: all affected is very likely because two mutations were",
    #   "detected at different reads and the tcn is 2 but", 
    #   "we are not sure if tumor is subclonal and therefore only ", 
    #   "wt-copies left")
    # log_obj[['tcn']] <- most_important_comb$tcn
    # log_obj[['aff_cp']] <- 
    #   most_important_comb$tcn-most_important_comb$wt_cp %>% 
    #   as.numeric()
    # log_obj[['class']] <- paste("comb:", most_important_comb$comb)
  } else if(nrow(all_comb)>3){
    vm("large number of variants", verbose)
    most_important_comb <- all_comb %>% 
      filter(wt_cp==min(all_comb$wt_cp)) %>%
      .[1,]
    log_obj[['score']] <- 1
    log_obj[['info']] <- 
      paste('unclear: very high number of mutations...', 
            'probably wt copies left but please check variant table')
    log_obj[['tcn']] <- most_important_comb$tcn
    log_obj[['aff_cp']] <- 
      most_important_comb$tcn-most_important_comb$wt_cp %>% 
      as.numeric()
    log_obj[['class']] <- paste("comb:", most_important_comb$comb)
  } else {
    vm("other", verbose)
    #print(names(all_comb))
    #print(all_comb %>% select(wt_cp))
    most_important_comb <- all_comb %>% 
      filter(wt_cp==min(all_comb$wt_cp)) %>%
      .[1,]
    log_obj[['score']] <- 1
    log_obj[['info']] <- 'wt-copies left' 
    log_obj[['tcn']] <- most_important_comb$tcn
    log_obj[['aff_cp']] <- 
      most_important_comb$tcn-most_important_comb$wt_cp %>% 
      as.numeric()
    log_obj[['class']] <- paste("comb:", 
                                most_important_comb$comb)            
  }
  vm("finished evaluation, returning results", verbose)
  log_obj[['info']] <- paste("phasing:", log_obj[['info']])
  return(log_obj)
}
#' @keywords internal
#' desicion tree if one variant already affects all copies in pre evaluation
eval_one_mut_affects_all <- function(df_gene, log_obj, printLog){
  concern_info <- df_gene %>% 
    filter(str_detect(pre_info, "all copies affected")) %>% .[1,]
  class <- get_classification(concern_info)
  origin <- ifelse(str_detect(concern_info$origin, "germ"),
                   "germ-",
                   "som-")
  log_obj[['score']] <- 2
  log_obj[['info']] <- concern_info$pre_info                 
  log_obj[['tcn']] <- concern_info$tcn
  log_obj[['aff_cp']] <- as.numeric(concern_info$aff_cp) 
  log_obj[['class']] <- paste0(origin, class)
  catt(printLog, 1, c("all copies affected by variant:", 
                      concern_info$mut_id))
  catt(printLog, 1, concern_info$pre_info)
  return(log_obj)
}
#' @keywords internal
#' if no variant affecta all copies
eval_one_mut <- function(df_gene, log_obj, printLog){
  log_obj[['score']] <- 1
  log_obj[['info']] <- df_gene$pre_info
  log_obj[['tcn']] <- df_gene$tcn
  log_obj[['aff_cp']] <- as.numeric(df_gene$aff_cp)
  log_obj[['class']] <- paste0(get_classification(df_gene))
  catt(printLog, 1, "1 variant detected:")
  catt(printLog, 1, df_gene$pre_info)
  return(log_obj)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_replace_all str_count str_detect
#' @importFrom dplyr as_tibble group_by relocate bind_rows left_join select tally mutate filter nth relocate
predict_zygosity_genewise <- function(GENE, df_all_mutations, bamDna, 
                                      bamRna, showReadDetail,
                                      printLog, purity, vcf, haploBlocks,
                                      phasedVcf, distCutOff, 
                                      verbose, logDir){
  if(printLog==TRUE){
    message(GENE)
  }
  gene <- pre_info <-  chr <-  pos <- wt_cp <-  mut_id <- status <- . <- 
    score <- comb <- dist <- tcn <- info <- n <- NULL 
  start_gene_eval <- Sys.time()
  pre_df_gene <- df_all_mutations %>% filter(gene==GENE)    
  df_gene_raw <- pre_df_gene%>%
    # all germline variants which are lost in the tumor are excluded
    filter(!str_detect(pre_info, 'variant lost in tumor'))
  ## check if variants at the same position are present
  n_vars_per_pos <- df_gene_raw %>%
    group_by(chr, pos) %>%
    tally() %>%
    filter(n>1)
  if(nrow(n_vars_per_pos)!=0){
    var_to_remove <- apply(n_vars_per_pos, 1, function(VAR){
      selected <- df_gene_raw %>% 
        filter(chr==VAR[["chr"]],
               pos==VAR[["pos"]]) %>%
        filter(wt_cp==min(.$wt_cp)) %>%
        ## if both vars have the same impac we just select the first
        .[1,]
      warning(
          "more than one variant detected at position:", 
          VAR[["chr"]], ":", VAR[["pos"]],
          "\n  maybe an indel with shifted annotation of position?",
          "\nSelecting variant of highest relevance by lowest left wt-copies:",
          selected$mut_id)
      to_remove <- df_gene_raw %>% 
        filter(chr==VAR[["chr"]],
               pos==VAR[["pos"]]) %>%
        filter(mut_id!=selected$mut_id) %>%
        select(chr, pos, mut_id) %>%
        mutate(remove=TRUE)
      return(to_remove)
    }) %>%
      bind_rows()
    df_gene <- df_gene_raw %>%
      left_join(var_to_remove, 
                by=c("chr"="chr", "pos"="pos", "mut_id"="mut_id")) %>%
      filter(is.na(remove)) %>%
      select(-remove)
  } else {
    df_gene <- df_gene_raw
  }
  if(!nrow(df_gene)==0){
    log_obj <- list(gene=GENE, n_mut=nrow(df_gene))
    all_comb <- NULL
    all_comb_to_export <- NULL
    #read_details <- NULL
    #snp_phasing <- NULL
    ## (1): check if any of the mutations already affects all alleles           
    if(str_detect(paste(df_gene$pre_info, collapse = ' '), 
                  'all copies affected')){
      eval_for_gene <- eval_one_mut_affects_all(df_gene, log_obj, printLog)
    ## (2): if only one variant is present  
    } else if(nrow(df_gene)==1){
      eval_for_gene <- eval_one_mut(df_gene, log_obj, printLog)
    ## (3): more than one variant present in gene   
    } else {  
      
      catt(printLog, 1, c(nrow(df_gene) ,
            "variants detected: Initializing haplotype phasing"))
      
      ## returns two tibbles: first one direct phasing combinations
      ## second one indirect phasing combinations
      full_phasing_result <- phase(df_gene, bamDna, bamRna, showReadDetail,
                        purity, vcf, haploBlocks, phasedVcf,
                        distCutOff, printLog, verbose, logDir)
      vm("full phasing done", verbose)
      all_comb <- full_phasing_result[[1]] %>%
        mutate(gene=GENE)
      
      #print(all_comb %>% select(gene, status, left_wt_copies))
      #print(names(all_comb))
      
      all_comb_to_export <- all_comb %>% 
        select(-score) %>%
        relocate(gene, comb, dist, tcn, 
                        status, left_wt_cp=wt_cp)
      
      eval_for_gene <- eval_phasing(all_comb, df_gene, log_obj, 
                                    printLog, verbose)
      
      #print(all_comb_to_export)
      # snp_phasing <- lapply(all_comb_raw, nth, n=3) %>%   
      #   bind_rows() %>%
      #   mutate(gene=GENE)
      # if(showReadDetail==TRUE){
      #   read_details <- lapply(all_comb_raw, nth, n=2) %>%   
      #     bind_rows() %>%
      #     mutate(gene=GENE)        
      # }
      
    } 
    vm("evaluation for gene done", verbose)
    #print(log_obj)
    df_reduced <- eval_for_gene %>% as_tibble() %>%
      mutate(status=ifelse(score==2,
                           "all_copies_affected",
                           "wt_copies_left"),
             info=str_replace_all(info, 
                                  " -> all copies affected| -> wt copies left",
                                  "")) 
    zygosity_gene <- list(pre_df_gene, 
                          list(df_reduced, 
                               all_comb_to_export)#, 
                          #read_details#, 
                          #snp_phasing
                          )
  } else {
    zygosity_gene <- list(pre_df_gene, NULL)
  } 
  if(printLog==TRUE){
    end_gene_eval <- Sys.time()
    message("~", round(
                    as.numeric(
                      difftime(
                        end_gene_eval, 
                        start_gene_eval,
                        units="secs")
                      ),3), "s") 
  }
  vm("returning result", verbose)
  return(zygosity_gene)
}
#' @keywords internal
#' description follows 
#' @importFrom dplyr case_when
get_classification <- function(data){
  case_when(
    data$class=="homdel" ~ "homdel",
    nchar(data$alt)==nchar(data$ref) ~ "snv",
    nchar(data$alt)>nchar(data$ref) ~ "ins",
    nchar(data$alt)<nchar(data$ref) ~ "del",
  ) %>%
    return()
}
#' @keywords internal
#' description follows
#' @importFrom dplyr case_when
define_class <- function(ref, alt){
  case_when(
    nchar(alt)==nchar(ref) ~ "snv",
    nchar(alt)>nchar(ref) ~ "ins",
    nchar(alt)<nchar(ref) ~ "del",
  ) %>%
    return()
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr as_tibble filter select mutate mutate_all
combine_uncovered_input_variants <- function(somSmallVars, germSmallVars,
                                             som_covered, germ_covered,
                                             templateGenes){
  mid <- seqnames <- start <- ref <- alt <- gene <- uncovered_som <- 
    uncovered_germ <- NULL
  if(!is.null(somSmallVars)){
    uncovered_som <-  as_tibble(somSmallVars) %>%
      filter(!mid %in% som_covered) %>%
      select(chr=seqnames, pos=start,ref, alt, gene) %>%
      mutate(origin="somatic") %>%
      mutate_all(.funs = as.character)
  } 
  if(!is.null(germSmallVars)){
    uncovered_germ <- as_tibble(germSmallVars) %>%
      filter(!mid %in% germ_covered) %>%
      select(chr=seqnames, pos=start,ref, alt, gene) %>%
      mutate(origin="germline") %>%
      mutate_all(.funs = as.character)
  }
  
  if(nrow(bind_rows(uncovered_som, uncovered_germ))==0){
    combined_uncovered <- NULL
  } else {
    combined_uncovered <- bind_rows(uncovered_som, uncovered_germ) %>%
      filter(gene %in% templateGenes)
  }
  return(combined_uncovered)
}
#' @keywords internal
check_opt_assgap <- function(assumeSomCnaGaps, ploidy){
  if(assumeSomCnaGaps==TRUE&is.null(ploidy)){
    warning("somatic CNA gaps can only be assumed if input ploidy is",
            "provided. Provide ploidy=2 to assume diploid case")
    return(FALSE)
  } else {
    return(assumeSomCnaGaps)
  }
}
#' @keywords internal
check_opt_incdel <- function(includeIncompleteDel, ploidy){
  if(includeIncompleteDel==TRUE&is.null(ploidy)){
    warning("Large scale deletions cannot be included without ploidy",
                "Please provide input ploidy. Provide ploidy=2",
                "to assume diploid case")
    return(FALSE)
  } else {
    return(includeIncompleteDel)
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr as_tibble group_by mutate ungroup filter
#' @importFrom purrr compact
combine_main_variant_tables <- function(df_germ, df_som, df_homdels,
                                        templateGenes, purity){
  gene <- . <- NULL
  df_all_mutations_unfiltered <- list(df_germ, df_som, df_homdels) %>% 
    compact() %>%
    Reduce(function(x,y)rbind(x,y),.) %>% 
    as_tibble() %>% 
    group_by(gene) %>% 
    mutate(mut_id=paste0("m", seq_len(length(gene))))%>% 
    ungroup() %>%
    mutate(purity=purity) 
  df_all_mutations <- df_all_mutations_unfiltered %>%
    filter(gene %in% templateGenes)
  if(nrow(df_all_mutations_unfiltered)!=nrow(df_all_mutations)){
    warning(
      abs(nrow(df_all_mutations_unfiltered)-nrow(df_all_mutations)),
      "variants in genes not found in reference genome annotation",
      "(geneModel). They are removed from the analysis.",
      "Please provide a annotation that contains every gene from the",
      "variant inputs"
    )
  }
  return(df_all_mutations)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr bind_rows group_by mutate relocate select ungroup relocate
bind_incdel_to_pre_eval <- function(df_incompletedels, df_all_mutations){
  purity <- gene <- mut_id <- NULL
  if(is.null(df_all_mutations)&is.null(df_incompletedels)){
    full_df_all_mutations <- NULL
  } else {
    if(is.null(df_incompletedels)){
      full_df_all_mutations <- df_all_mutations %>%
        select(-purity) %>%
        relocate(gene, mut_id, 2:13)
    } else {
      full_df_all_mutations <- bind_rows(
        df_all_mutations,
        df_incompletedels %>% 
          unique() %>% 
          mutate(mut_id=NA)
      ) %>% group_by(gene) %>% 
        mutate(mut_id=paste0("m", seq_len(length(gene))))%>% 
        ungroup()%>%
        select(-purity) %>%
        relocate(gene, mut_id, 2:13)       
    }    
  }

  return(full_df_all_mutations)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows mutate select
bind_incdel_to_final_eval <- function(df_incompletedels, final_output){
  gene <- 0
  if(is.null(final_output)&is.null(df_incompletedels)){
    full_output <- NULL
  } else {
    if(is.null(df_incompletedels)){
      full_output <- final_output
    } else {
      full_output <- bind_rows(
        final_output,
        df_incompletedels %>%
          filter(!gene %in% final_output$gene) %>%
          select(gene) %>%
          unique() %>%
          mutate(status="wt_copies_left",
                 info="somatic-incompletedel")
      )       
    }    
  }
  return(full_output)
}
#' @importFrom IRanges subsetByOverlaps
#' @export
predict_per_variant <- function(purity, 
                             sex,
                             somCna, 
                             geneModel=NULL,
                             somSmallVars=NULL, 
                             germSmallVars=NULL,
                             ploidy=NULL, 
                             colnameTcn=NULL,
                             colnameCnaType=NULL,
                             includeHomoDel=TRUE,
                             includeIncompleteDel=TRUE,
                             assumeSomCnaGaps=FALSE,
                             byTcn=TRUE,
                             is_pre_eval=TRUE
){
  status <- info <- wt_cp <- . <- df_homdels <- df_all_mutations <- gene <-
    final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- 
    combined_snp_phasing <- NULL
  ## check input for valid format
  somSmallVars <- check_gr_small_vars(somSmallVars, "somatic")
  germSmallVars <- check_gr_small_vars(germSmallVars, "germline")
  purity <- check_purity(purity)
  ploidy <- check_ploidy(ploidy)
  sex <- check_sex(sex)
  geneModel <- check_gr_gene_model(geneModel, is_pre_eval)
  assumeSomCnaGaps <- check_opt_assgap(assumeSomCnaGaps, ploidy)
  somCna <- check_somCna(somCna, geneModel, sex, ploidy, 
                         assumeSomCnaGaps, colnameTcn, 
                         colnameCnaType) 
  includeIncompleteDel <- check_opt_incdel(includeIncompleteDel, ploidy)
  includeHomoDel <- check_opt_incdel(includeHomoDel, ploidy)
  if(is.null(geneModel)){
    ## only small variants can be evaulated
    if(includeIncompleteDel==TRUE|includeHomoDel==TRUE){
      includeIncompleteDel <- FALSE
      includeHomoDel <- FALSE
      warning("To include large deletions the geneModel muist be provided.",
              "IncludeHomoDel and IncludeIncompleteDel will be FALSE")
    }
    templateGenes <- c(somSmallVars$gene, germSmallVars$gene)
  } else {
    templateGenes <- geneModel$gene
  }
  
  ## get variants not covered by CNV input
  if(!is.null(somSmallVars)){
    gr_som_cov <- subsetByOverlaps(somSmallVars, somCna) 
    som_covered <- gr_som_cov$mid
  }
  if(!is.null(germSmallVars)){
    gr_germ_cov <- subsetByOverlaps(germSmallVars, somCna) 
    germ_covered <- gr_germ_cov$mid
  }
  combined_uncovered <- combine_uncovered_input_variants(somSmallVars, 
                                                         germSmallVars,
                                                         som_covered,
                                                         germ_covered,
                                                         templateGenes)
  ## preapre input data for prediction (pre-evaluation)
  df_germ <- prepare_germline_variants(
    gr_germ_cov, 
    somCna, purity, sex)
  df_som <- prepare_somatic_variant_table(
    gr_som_cov, 
    templateGenes, 
    somCna, purity, sex)
  df_homdels <- extract_all_dels_of_sample(somCna, geneModel, 
                                           "homdel", byTcn, sex, ploidy,
                                           includeHomoDel)
  df_incompletedels <- extract_all_dels_of_sample(somCna, geneModel, 
                                                  "incompletedel", TRUE, sex, 
                                                  ploidy, includeIncompleteDel)
  ## predict zygosity for each gene
  if(!is.null(df_germ)|!is.null(df_som)|!is.null(df_homdels)){
    df_all_mutations <- combine_main_variant_tables(df_germ, df_som, 
                                                    df_homdels,
                                                    templateGenes, purity)
  }
  if(is_pre_eval==TRUE){
    return(list(evaluation_per_variant= 
                  bind_incdel_to_pre_eval(df_incompletedels, df_all_mutations),
                combined_uncovered=combined_uncovered
      
    ))
  } else {
    return(list(evaluation_per_variant=df_all_mutations,
                df_incompletedels=df_incompletedels,
                combined_uncovered=combined_uncovered))    
  }

}
vm <- function(mes, verbose=FALSE){
  if(verbose==TRUE){
    #message(paste0("\n",mes))
    print(mes)
  }
}
#' predicts zygosity of a set of genes of a sample
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param ploidy ploidy of the sample (numeric value)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param somCna GRanges object containing all genomic regions with annotated 
#' total copynumber and cna_type as metadata columns. The total-copynumber 
#' column should be named "tcn" but also some other commonly used names. 
#' It should contain numeric values or characters that can be converted to 
#' numeric values. The cna_type column must contain the information about 
#' loss of heterozygosity (LOH). Therefore the term "LOH" must be explicitely 
#' mentioned in the column. If a genomic region is not present in the object, 
#' it will be taken as heterozygous with neutral TCN of 2. 
#' @param somSmallVars GRanges object containing all somatic small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), 
#' alternative base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF). 
#' If the object is not provided the tool assumes there are no somatic small 
#' variants.
#' @param germSmallVars GRanges object containing all germline small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), alternative 
#' base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF)
#' If the object is not provided the tool assumes there are no germline small 
#' variants.
#' @param geneModel GRanges object containing the gene-annoattion of 
#' the used reference genome with metadata column of the gene name (gene)
#' @param bamDna path to bam-file
#' @param bamRna optional; path to rna file (bam format)
#' @param includeHomoDel default = TRUE; if FALSE homozygous deleteions are 
#' excluded
#' @param includeIncompleteDel default = TRUE; if FALSE heterzygous deleteions 
#' are excluded
#' @param showReadDetail default = FALSE; if TRUE a table is added to the 
#' output,
#' @param printLog default = FALSE; if TRUE the gene which is evaluated is 
#' printed in console, 
#' containing the query-name of each read which was used to perform 
#' haplotype-phasing and the info into which class it was assigned.
#' @param assumeSomCnaGaps (logical, default=FALSE) Only required if the somCna
#' object lacks copy number information for genomic segments on which small 
#' variants are detected. By default, variants in such regions will be excluded 
#' from the analysis as required information about the copy number is missing. 
#' These variants will be attached to the final output list in a separate 
#' tibble. To include them, this flag must be set TRUE and the ground ploidy 
#' must be given as an input. This ground ploidy will then be taken as tcn in 
#' the missing regions. If no ploidy is given the tool will assume the ground 
#' ploidy of 2 when this flag is TRUE.
#' @param byTcn logical, default=TRUE; optional if includeHomoDel or 
#' includeIncompleteDelS is TRUE. If FALSE the tool will not use tcn as a 
#' criterion to assign large deletions. It will use the cna_type column and 
#' check for indicating strings like HOMDEL/HomoDel/DEL. Some commonly used 
#' strings are covered. It is recommended to leave this flag TRUE
#' @param colnameTcn character indicating the name of the metadata containing 
#' the tcn information in the somCna object. If not provided the tool tries to 
#' detect the column according to default names
#' @param colnameCnaType character indicating the name of the metadata 
#' containing cna type information in the somCna object. 
#' If not provided the tool tries to detect the column according to default 
#' names
#' @param vcf character; path to variant call file (.vcf.gz format). 
#' Will be used (if provided)
#' for extended SNP phasing if variants on the same gene are too far away from
#' each other for direct haplotype phasing
#' @param distCutOff numeric, default=5000; if input vcf is provided and SNP
#' phasing is performed, this will limt the distance at which the SNP phasing
#' should not be tried anymore. As the probability of finding overlapping reads
#' at such a long distance is very low and the runtime will increase
#' exponentially.
#' @return A list of dataframes. Those are the evaluation per variant, 
#' the evaluation per gene and, if performed, the info about the 
#' haplotype-phasing.
#' @examples
#' cnvs  = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr = "chr17",
#'     start = c(170060, 34520990),
#'     end = c(34520990, 83198614),
#'     tcn = c(2, 1),
#'     cna_type = c("neutral", "LOH")
#'   )
#' )
#' somatic_vars = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr="chr17",
#'     start = 7675088,
#'     end = 7675088,
#'     ref = "C",
#'     alt = "T",
#'     af = 0.65,
#'     gene = "TP53" 
#'   )
#' )
#' germline_vars = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr="chr17",
#'     start = 41771694,
#'     end = 41771694,
#'     ref = "GTGT",
#'     alt = "G",
#'     af = 0.95,
#'     gene = "JUP" 
#'   )
#' )
#' reference = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr = "chr17",
#'     start = c(7661778, 41754603),
#'     end = c(7687538, 41786931),
#'     gene = c("TP53", "JUP")
#'   )
#' )
#' sex = "female"
#' purity = 0.9
#' bamfile <- system.file("extdata", "ZP_example.bam", 
#'   package = "ZygosityPredictor")
#' predict_zygosity(purity = purity, sex = sex, 
#'   somCna = cnvs,
#'   somSmallVars = somatic_vars,
#'   germSmallVars = germline_vars,
#'   geneModel = reference,
#'   bamDna = bamfile
#' )
#' @importFrom stringr %>%
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr compact
#' @importFrom dplyr bind_rows nth select tibble
#' @export
predict_zygosity <- function(purity, 
                             sex,
                             somCna, 
                             geneModel,
                             bamDna,
                             somSmallVars=NULL, 
                             germSmallVars=NULL,
                             bamRna=NULL,
                             ploidy=NULL, 
                             colnameTcn=NULL,
                             colnameCnaType=NULL,
                             includeHomoDel=TRUE,
                             includeIncompleteDel=TRUE,
                             showReadDetail=FALSE,
                             printLog=FALSE,
                             assumeSomCnaGaps=FALSE,
                             byTcn=TRUE,
                             vcf=NULL,
                             haploBlocks=NULL,
                             phasedVcf=NULL,
                             distCutOff=5000,
                             verbose=FALSE,
                             logDir=NULL
                             ){
  status <- info <- wt_cp <- . <- df_homdels <- df_all_mutations <- gene <-
    final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- 
    combined_snp_phasing <- evaluation_per_gene <- NULL
  #vm("initializing evaluation per variant", verbose)
  evaluation_per_variant_pre <- predict_per_variant(purity, 
                                                sex,
                                                somCna, 
                                                geneModel,
                                                somSmallVars, 
                                                germSmallVars,
                                                ploidy, 
                                                colnameTcn,
                                                colnameCnaType,
                                                includeHomoDel,
                                                includeIncompleteDel,
                                                assumeSomCnaGaps,
                                                byTcn,
                                                is_pre_eval=FALSE)
  evaluation_per_variant <- evaluation_per_variant_pre$evaluation_per_variant
  if(!is.null(evaluation_per_variant)){
    df_all_mutations <- evaluation_per_variant
    if(!nrow(df_all_mutations)==0){
      bamDna <- check_bam(bamDna)
      bamRna <- check_rna(bamRna)
      logDir <- check_logDir(logDir)
      vcf <- check_vcf(vcf) 
      haploBlocks <- check_haploblocks(haploBlocks)
      #vm("initializing evaluation per gene", verbose)
      result_list <- lapply(
        unique(df_all_mutations$gene), 
        predict_zygosity_genewise, 
        df_all_mutations, 
        bamDna,
        bamRna,
        showReadDetail,
        printLog,
        purity,
        vcf,
        haploBlocks,
        phasedVcf,
        distCutOff, 
        verbose,
        logDir)
      #print(result_list)
      pre_scoring <- lapply(result_list, nth, n=1) %>% 
        bind_rows()
      final_result_list <- lapply(result_list, nth, n=2) %>% compact()
      if(length(final_result_list)!=0){
        final_phasing_info <- lapply(final_result_list, nth, n=2) %>% 
          compact() %>% 
          bind_rows() 
        final_output <- lapply(final_result_list, nth, n=1) %>% 
          bind_rows() %>%
          select(gene, status, info)
        # combined_read_details <- lapply(result_list, nth, n=3) %>% compact() %>%
        #   bind_rows()
        # combined_snp_phasing <- lapply(result_list, nth, n=4) %>% compact() %>%
        #   bind_rows()
      } 
    }
    evaluation_per_gene <- bind_incdel_to_final_eval(
      evaluation_per_variant_pre$df_incompletedels,
      final_output)
  }
  result_list <- list(
    eval_per_variant=evaluation_per_variant,
    eval_per_gene=evaluation_per_gene,
    phasing_info=final_phasing_info,
    #readpair_info=combined_read_details,
    uncovered_input=evaluation_per_variant_pre$combined_uncovered#,
    #ext_snp_phasing=combined_snp_phasing
  ) %>%
    compact()
  return(result_list)
} 