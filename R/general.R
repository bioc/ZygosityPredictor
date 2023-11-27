func_end <- function(verbose=NULL){
  mes <- "end"
  if(!is.null(verbose)){
    vm(mes, verbose, -1)
  } else {
    vm(mes, vb, -1)
  }
}
func_start <- function(verbose=NULL){
  #vm(as.character(sys[1]), verbose, 1)
  mes <- paste0(as.character(sys.call(-1)[1]))
  if(!is.null(verbose)){
    vm(mes, 
       verbose, 1)
  } else {
    vm(mes, vb, 1)
  }
  
}
increment_log_depth <- function() {
  log_depth <<- log_depth + 1
}
decrement_log_depth <- function() {
  log_depth <<- log_depth - 1
}
increment_call_depth <- function() {
  call_depth <<- call_depth + 1
}
decrement_call_depth <- function() {
  call_depth <<- call_depth - 1
}
store_log <- function(geneDir, obj, file){
  if(!is.null(geneDir)){
    if(!is.null(obj)){
      write_tsv(obj, 
                file=file.path(geneDir, 
                               file))      
    }
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
  func_start(verbose)
  qname.first <- . <- NULL
  ## importFrom dplyr tibble filter
  ref_pos1 <- as.numeric(pos1)
  ref_pos2 <- as.numeric(pos2)
  ref_chr1 <- as.character(chr1)
  ref_chr2 <- as.character(chr2)
  #print(pos1)
  #print(pos2)
  if(pos1>pos2){
    #stop("avoid changing position order")
    ## exchange does not matter here
    ref_pos1 <- as.numeric(pos2)
    ref_pos2 <- as.numeric(pos1)
    ref_chr1 <- as.character(chr2)
    ref_chr2 <- as.character(chr1)
  }
 #print(ref_chr1)
 #print(ref_chr2)
  ref_gr1 <- GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1)
  ref_gr2 <- GRanges(seqnames = ref_chr2, 
                                   ranges = ref_pos2)
  ## now load all reads/read-pairs that cover the position of the first variant
  #vm("loading reads", verbose, 1)
  all_covering_read_pairs <- readGAlignmentPairs(
    bamDna,
    param=ScanBamParam(
      which=GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1),
      what=c("qname","seq", "cigar", "mapq", "qual")
    )) 
  if(length(all_covering_read_pairs)==0){
    filtered_reads <- tibble()
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
    
  }
  func_end(verbose)
  return(filtered_reads)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom purrr compact
#' @importFrom dplyr summarize pull
insert_missing_cnv_regions <- function(somCna, geneModel, sex, ploidy){
  vm(as.character(sys.call()[1]), verbose, 1)
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
  func_end(verbose)
  return(comb_somCna)
}

#' calc_left_wt_copies <- function(mtcn, nstatus, aff_copies1, aff_copies2){
#'   left_wt_copies <- ifelse(nstatus==2,
#'                            mtcn-sum(as.numeric(aff_copies1), 
#'                                     as.numeric(aff_copies2)),
#'                            mtcn-max(as.numeric(aff_copies1), 
#'                                     as.numeric(aff_copies2)))
#'   return(as.numeric(left_wt_copies))
#' }
#' 
#' 
#' aggregate_probs <- function(ps){
#'   if(length(ps)>1){
#'     np <- ps[1]
#'     for (n in seq(1, length(ps)-1, 1)){
#'       np <- np+(1-np)*ps[n+1]
#'     }
#'     return(np)    
#'   } else {
#'     return(ps) 
#'   }
#' }
#' 
#' extract_subseq <- function(ref_pos, element_mut, parsed_read, 
#'                            length_indel, string){
#'   pos_in_seq <- ref_pos-element_mut$map_start+1
#'   subseq_in_element <- str_sub(element_mut[[string]], start=pos_in_seq, 
#'                                end=-1)
#'   subseq_in_following_elements <- 
#'     parsed_read[
#'       which(as.numeric(parsed_read$id)>as.numeric(element_mut$id)),] %>%
#'     pull(all_of(string)) %>%
#'     na.omit() %>%
#'     paste(collapse="")
#'   full_remaining <- paste0(subseq_in_element, subseq_in_following_elements)
#'   return(str_sub(full_remaining, start=1, end=length_indel))
#' }
#' 
#' 
#' extract_snv <- function(ref_pos, element_mut){
#'   ## get cigar element in which reference position in located
#'   pos_in_seq <- ref_pos-element_mut$map_start+1
#'     base <- 
#'       str_sub(element_mut$seq,
#'               start=pos_in_seq,
#'               end=pos_in_seq) 
#'     qual <- 
#'       str_sub(element_mut$qual,
#'               start=pos_in_seq,
#'               end=pos_in_seq)  
#'   return(tibble(base=base, qual=qual, info="mapped", len_indel=NA,
#'                 exp_indel=NA))
#' }
#' 
#' extract_insertion <- function(ref_pos, parsed_read, element_mut, length_ins){
#'   ## insertion are always indicated by the I type. The reference position tells
#'   ## the position in front of the I segment, which means we are looking for
#'   ## the subsequent one type
#'   expected_indel_detected <- FALSE
#'   len <- NA
#'   if(as.numeric(element_mut$id)==nrow(parsed_read)){
#'     if(ref_pos==element_mut$map_end){
#'       info <- 
#'         "no insertion detected: ref pos is last mapped base of read" 
#'     } else {
#'       info <- 
#'         "no insertion detected: ref pos mapped in last element of parsed read"      
#'     }
#'     base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
#'     #qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
#'   } else if("I" %in% parsed_read$type){
#'     ## insertion detected... check if it is the matching one
#'     if(ref_pos==element_mut$map_end){
#'       ## if insertion is at refrence position, the subsequent position is an a
#'       ## new elemnt labelled as I
#'       next_element <- parsed_read[as.numeric(element_mut$id)+1,]
#'       if(next_element$type=="I"){
#'         ## insertion detected at position.. extract inserted sequence
#'         base <- paste0(
#'           str_sub(element_mut$seq, start=-1, end=-1),
#'           next_element$seq
#'         )
#'       #  qual <- paste0(
#'        #   str_sub(element_mut$qual, start=-1, end=-1),
#'         #  next_element$qual
#'         #)
#'         info <- "insertion detected"
#'         len <- next_element$width
#'         expected_indel_detected <- TRUE
#'       } else {
#'         info <- paste(next_element$type, "detected - wrong class")
#'         base <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
#'                                "seq")
#'         #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
#'          #                      "qual")
#'       }
#'     } else {
#'       info <- "no insertion detected: ref pos not mapped to end of element"
#'       base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
#'      # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
#'     }
#'   } else {
#'     ## no inserttion detected.. extract base at position to get base quality
#'     info <- "no insertion detected: no I in cigar"
#'     base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
#'    # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
#'   }
#'   ## qual is set to null, as we only take the mapping quality for indels
#'   return(tibble(base=base, qual=0, info=info, len_indel=len, 
#'                 exp_indel=expected_indel_detected))
#' }
#' 
#' extract_deletion <- function(ref_pos, parsed_read, element_mut, length_del){
#'   expected_indel_detected <- FALSE
#'   len <- NA
#'   if(as.numeric(element_mut$id)==nrow(parsed_read)){
#'     if(ref_pos==element_mut$map_end){
#'       info <- 
#'         "no deletion detected: ref pos is last mapped base of read" 
#'     } else {
#'       info <- 
#'         "no delection detected: ref pos mapped in last element of parsed read"      
#'     }
#'     base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
#'    # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
#'   } else if("D" %in% parsed_read$type){
#'     ## insertion detected... check if it is the matching one
#'     if(ref_pos==element_mut$map_end){
#'       ## if insertion is at refrence position, the subsequent position is an a
#'       ## new elemnt labelled as D
#'       next_element <- parsed_read[as.numeric(element_mut$id)+1,]
#'       if(next_element$type=="D"){
#'         expected_indel_detected <- TRUE
#'         info <- "deletion detected"
#'         len <- next_element$width
#'       } else {
#'         info <- "no deletion detected: subsequent element not I"
#'       }
#'     } else {
#'       info <- "no insertion detected: ref pos not mapped to end of element"
#'     }      
#'     base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
#'                            "seq")
#'     #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
#'      #                      "qual")
#'   } else {
#'     ## no inserttion detected.. extract base at position to get base quality
#'     info <- "no insertion detected: no I in cigar"
#'     base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
#'                            "seq")
#'     #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
#'      #                      "qual")
#'   }
#'   ## qual is set to null, as we only take the mapping quality for indels
#'   return(tibble(base=base, qual=0, info=info, len_indel=len,
#'                 exp_indel=expected_indel_detected))
#' }
#' 
#' cigar_element <- function(parsed_read, ref_pos){
#'   #func_start()
#'   cigel <- parsed_read[which(between(rep(ref_pos, nrow(parsed_read)), 
#'                                      parsed_read$map_start, 
#'                                      parsed_read$map_end)),] %>%
#'     .[1,]
#'   #func_end()
#'   return(cigel)
#' }
#' 
#' extract_base_at_refpos <- function(parsed_read, ref_pos, class, ref_alt, 
#'                                    ref_ref){
#'   #func_start()
#'   ## extract elemnt containing reference position from parsed read
#'   element_mut <- cigar_element(parsed_read, ref_pos)
#'   ## if position is inside N or D type, the read does not cover the position
#'   ## N means for example skipped exon in RNA
#'   ## D wouldmean the position is deleted
#'   ## real deletions are still mapped in the M part
#'   if(element_mut$type %in% c("N", "D")){
#'     base_info <- tibble(base=NA, qual=NA, info="position skipped/deleted", 
#'                         len_indel=NA, exp_indel=NA)
#'   } else if(class %in% c("snv", "snp")){
#'     base_info <- extract_snv(ref_pos, element_mut)
#'   } else if(class=="ins"){
#'     base_info <- extract_insertion(ref_pos, parsed_read, element_mut, 
#'                                    nchar(ref_alt)-nchar(ref_ref))
#'   } else if(class=="del"){
#'     base_info <- extract_deletion(ref_pos, parsed_read, element_mut, 
#'                                   1+abs(nchar(ref_ref)-nchar(ref_alt)))
#'   } else {
#'     stop("class of variant needs to be provided")
#'   }
#'   #func_end()
#'   return(base_info %>% mutate(class=class, mapq=element_mut$mapq))
#' }
#' 
#' #' @keywords internal
#' #' @importFrom stringr %>% str_replace str_sub
#' #' @importFrom dplyr between case_when rowwise filter
#' check_mut_presence <- function(read_structure, ref_pos,
#'                                   ref_alt, ref_ref, ref_class, verbose=FALSE){
#'   map_start <- map_end <- . <- NULL
#'   ## check if read is spanned out:
#'   element_mut <- read_structure %>% rowwise() %>%
#'       filter(between(ref_pos, map_start, map_end)) %>%
#'       .[1,]
#'   if(element_mut$type %in% c("N", "D")){
#'     mut_at_read <- -10
#'   } else {
#'     if(ref_class %in% c("snv", "snp")){
#'       pos_in_seq <- ref_pos-element_mut$map_start+1
#'       base <- 
#'         str_sub(element_mut$seq,
#'                 start=pos_in_seq,
#'                 end=pos_in_seq) 
#'       qual <- 
#'         str_sub(element_mut$qual,
#'                 start=pos_in_seq,
#'                 end=pos_in_seq)
#'       mut_at_read <- case_when(
#'         base==ref_alt ~ 1,
#'         base==ref_ref ~ 0,
#'         TRUE ~ -2
#'       ) 
#'     } else if(ref_class=="ins"){
#'       cat_at_pos <- 
#'         read_structure[which(read_structure$map_end==ref_pos)+1,] %>%
#'         .[which(.$type=="I"),] 
#'       if(nrow(cat_at_pos)==0){
#'         mut_at_read <- 0
#'         base <- NA
#'         qual <- NA
#'       } else {  
#'         length_ins <- str_replace(ref_alt, paste0("^", ref_ref), "") %>%
#'           nchar()
#'         if(length_ins==cat_at_pos[1,]$width){
#'           mut_at_read <- 1
#'         } else {
#'           mut_at_read <- -2
#'         }
#'         base <- cat_at_pos[1,]$seq
#'         qual <- cat_at_pos[1,]$qual
#'       }
#'     } else if(ref_class=="del"){
#'       cat_at_pos <- 
#'         read_structure[which(read_structure$map_end==ref_pos)+1,]%>%
#'         .[which(.$type=="D"),]
#'       
#'       if(nrow(cat_at_pos)==0){
#'         mut_at_read <- 0
#'         base <- NA
#'         qual <- NA
#'       } else {  
#'         length_del <- str_replace(ref_ref, paste0("^", ref_alt), "") %>%
#'           nchar()
#'         base <- NA
#'         qual <- NA
#'         if(length_del==cat_at_pos[1,]$width){
#'           mut_at_read <- 1
#'         } else {
#'           mut_at_read <- -2
#'         }
#'       }
#'     }
#'   }
#'   return(tibble(base=base,
#'                 qual=qual,
#'                 mut_at_read=mut_at_read))
#' }
#' parse_cigar <- function(bam, qname){
#'   #func_start()
#'   paired_reads <- bam[which(bam$qname==qname)] 
#'   read_start <- GenomicAlignments::start(paired_reads)
#'   seq <- as.character(paired_reads$seq)
#'   cigar <- paired_reads$cigar
#'   qual <- as.character(paired_reads$qual)
#'   mate <- c(1,2)
#'   ## parse cigar string according to query
#'   cigq <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar,
#'                                                         with.ops = T) 
#'   ## parse cigar string according to reference
#'   cigr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar,
#'                                                             with.ops = F) 
#'   raw_cigs <- lapply(mate, function(i){
#'     cigq_start <- GenomicAlignments::start(cigq[[i]])
#'     cigq_end <- GenomicAlignments::end(cigq[[i]])
#'     cigq_len <- length(cigq_start)
#'     cigr_start <- GenomicAlignments::start(cigr[[i]])
#'     cigr_end <- GenomicAlignments::end(cigr[[i]])
#'     seqq <- str_sub(seq[i],
#'                     start=cigq_start,
#'                     end=cigq_end) %>% na_if("")
#'     qualq <- str_sub(qual[i],
#'                      start=cigq_start,
#'                      end=cigq_end) %>% na_if("")
#'     width=pmax(GenomicAlignments::width(cigq[[i]]), 
#'                GenomicAlignments::width(cigr[[i]]))
#'     map_start=read_start[i]+cigr_end-width
#'     map_end=read_start[i]+cigr_end-1
#'     raw_cigs_new <- data.frame(
#'       width=width,
#'       type=names(cigq[[i]]),
#'       seq=seqq,
#'       qual=qualq,
#'       map_start=map_start,
#'       map_end=map_end,
#'       start=cigr_start,
#'       end=cigr_end,
#'       mate=mate[i],
#'       mapq=paired_reads$mapq[i],
#'       origin=paired_reads$origin[i]
#'     )
#'   }) %>%
#'     bind_rows() %>%
#'     mutate(id=c(1:nrow(.)))
#'   #func_end()
#'   return(raw_cigs)
#' }
#' evaluate_base <- function(base_info, ref_alt, ref_ref){
#'   #func_start()
#'   if(base_info$class %in% c("snv", "snp")){
#'     detected <- case_when(
#'       base_info$base==ref_alt ~ 1,
#'       base_info$base==ref_ref ~ 0,
#'       TRUE ~ -2
#'     )
#'   } else if(base_info$class=="ins"){
#'     detected <- case_when(
#'       base_info$base==ref_alt ~ 1,
#'       base_info$base==ref_ref ~ 0,
#'       TRUE ~ -2
#'     )
#'   } else {
#'     detected <- case_when(
#'       base_info$len_indel==abs(nchar(ref_alt)-nchar(ref_ref)) ~ 1,
#'       is.na(base_info$len_indel) ~ 0,
#'       TRUE ~ -2
#'     )
#'   }
#'   #func_end()
#'   return(detected)
#' }
#' #' @keywords internal
#' #' description follows
#' #' @importFrom stringr %>%
#' #' @importFrom tibble column_to_rownames rownames_to_column
#' #' @importFrom purrr set_names
#' #' @importFrom dplyr case_when bind_rows as_tibble mutate
#' #' @importFrom GenomicRanges seqnames
#' core_tool <- function(qname, bam,
#'                       ref_pos1, ref_pos2,
#'                       ref_alt1, ref_alt2,
#'                       ref_ref1, ref_ref2,
#'                       ref_class1, ref_class2,
#'                       verbose=FALSE, version="old"){
#'   #vm("core_tool", verbose, 1)
#'   . <- NULL
#'   ## parse read according to cigar string
#'   parsed_read <- parse_cigar(bam, qname)
#'   ## extract base at reference position
#'   base_info1 <- extract_base_at_refpos(parsed_read, ref_pos1, ref_class1, 
#'                                            ref_alt1, ref_ref1)
#'   base_info2 <- extract_base_at_refpos(parsed_read, ref_pos2, ref_class2, 
#'                                            ref_alt2, ref_ref2)
#'   ## assign final status to read
#'   if(is.na(base_info1$base)|is.na(base_info2$base)){
#'     final_assignment <- "skipped"
#'   } else {
#'     mut1_in_read <- evaluate_base(base_info1, ref_alt1, ref_ref1)
#'     mut2_in_read <- evaluate_base(base_info2, ref_alt2, ref_ref2)
#'     final_assignment <- case_when(
#'       sum(mut1_in_read, mut2_in_read)==2 ~ "both",
#'       sum(mut1_in_read, mut2_in_read)==0 ~ "none",
#'       mut1_in_read==1 ~ "mut1",
#'       mut2_in_read==1 ~ "mut2",
#'       TRUE ~ "dev_var"
#'     )    
#'   }
#'   #func_end()
#'   return(c(qname=qname, result=final_assignment,
#'            origin=unique(parsed_read$origin),
#'            baseq1=base_info1$qual,
#'            mapq1=base_info1$mapq,
#'            baseq2=base_info2$qual,
#'            mapq2=base_info2$mapq
#'            )
#'          )
#' }

#' @keywords internal
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate select mutate_at rowwise
#' @importFrom stringr str_detect
merge_sCNAs <- function(obj, somCna, verbose){
  func_start(verbose)
  alt <- af <- tcn <- tcn_assumed <- cna_type <- gene <- ref <- NULL
  merged <- obj %>%
    mergeByOverlaps(somCna) %>% as_tibble() %>%
    mutate(cna_type=ifelse(str_detect(cna_type, "LOH"), "LOH", "HZ")) %>%
    select(chr=1, pos=2, gene, ref, alt, af, tcn, cna_type, all_imb, 
           tcn_assumed) %>%
    mutate_at(.vars = c("af", "tcn"), .funs=as.numeric) %>%
    rowwise() #%>%
  func_end(verbose)
  return(merged)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate filter
prepare_somatic_variant_table <- function(somSmallVars, templateGenes, 
                                          somCna, purity, sex, verbose){
  func_start(verbose)
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- 
    aff_cp <- wt_cp <- . <- tcn_assumed <- NULL
  if(!is.null(somSmallVars)){
    reduced_to_templateGenes <- 
      somSmallVars[which(somSmallVars$gene %in% templateGenes)]
    if(length(reduced_to_templateGenes)>0){
      tbl_prepared_variants <- reduced_to_templateGenes %>%
        merge_sCNAs(., somCna, verbose) %>%
        
        ################################
        mutate(
          origin="somatic",
          class=define_class(ref, alt),
          aff_cp = aff_som_copies(chr, af, tcn, purity, sex),
          wt_cp = tcn-aff_cp,
          vn_status=case_when(
            (isTRUE(str_detect(cna_type,'LOH'))|chr=="X"&sex=="male")&
              wt_cp<=0.5 ~ 2,
            TRUE ~ 1
          ),
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
        )      
    } else {
      tbl_prepared_variants <- NULL
    }
  } else {
    tbl_prepared_variants <- NULL
  }
  func_end(verbose)
  return(tbl_prepared_variants)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr compact
#' @importFrom dplyr filter as_tibble select bind_rows mutate
extract_all_dels_of_sample <- function(somCna, geneModel, DEL_TYPE, 
                                       byTcn, sex, ploidy, include_dels, 
                                       verbose){
  func_start(verbose)
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
      #return(NULL)
      full_CNAs <- NULL
    } else {
      merged_CNAs <- apply(relevant_CNAs, 1, function(CNV){
        GR_single_CNV <- GRanges(as_tibble(t(CNV)))
        suppressWarnings(
        rel <- subsetByOverlaps(geneModel, GR_single_CNV) %>%
          as_tibble()
        )
        if(nrow(rel)==0){
          #return(NULL)
          full_CNAs <- NULL
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
                   all_imb=NA,
                   aff_cp=ploidy-as.numeric(CNV[["tcn"]]),
                   wt_cp=as.numeric(CNV[["tcn"]]),
                   pre_info=pre_info,
                   vn_status=ifelse(DEL_TYPE=="homdel", 2, 1)
            ) %>%
            return()
        }
      })%>% compact() %>% bind_rows() %>% unique()
      if(nrow(merged_CNAs)==0){
        full_CNAs <- NULL
        #return(NULL)
      } else {
        full_CNAs <- merged_CNAs
      }
    }
  } else {
    #return(NULL)
    full_CNAs <- NULL
  }
  func_end(verbose)
  return(full_CNAs)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom dplyr mutate select
prepare_germline_variants <- function(germSmallVars, somCna, purity, sex, 
                                      verbose){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- aff_cp <- 
    origin <- pos <- wt_cp <- pre_info <- . <- tcn_assumed <- NULL
  func_start(verbose)
  if(is.null(germSmallVars)){
    df_germ <- NULL
  } else {
    df_germ <- germSmallVars %>% 
      merge_sCNAs(., somCna, verbose) %>%
      mutate(
        origin="germline",
        class=define_class(ref, alt),
        aff_cp = aff_germ_copies(chr, af, tcn, purity, sex, 1),
        wt_cp = tcn-aff_cp,
        vn_status=case_when(
          isFALSE(cna_type=='LOH') ~ 1,
          af>=0.5 ~ 2,
          TRUE ~ -1 
        ),
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
             all_imb,
             aff_cp,
             wt_cp,
             pre_info,
             vn_status
      ) 
  }
  func_end(verbose)
  return(df_germ)
}


#' @keywords internal
#' description follows
check_for_overlapping_reads <- function(bamDna, bamRna,
                                        ref_chr1, 
                                        ref_chr2, 
                                        ref_pos1, 
                                        ref_pos2, 
                                        verbose=FALSE){
  func_start()
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
  func_end()
  return(c(dna_bam, rna_bam))
}
ascii_to_dec <- function(ascii_encoded){
  if(is.na(ascii_encoded)){
    return(NA)
  } else {
    exp <- (as.integer(charToRaw(ascii_encoded))-33)/(-10)
    if(length(exp)>1){
      exp <- mean(exp)
    }
    # Convert the QUAL string to raw hexadecimal
    return(10^exp)    
  }

}


print_tibble <- function(tbl_in){
  res <- paste(as.character(knitr::kable(tbl_in)), collapse="\n")
  return(res)
}
split_genomic_range <- function(gr, exclude_positions) {
  # Sort the exclude_positions and ensure they are unique
  exclude_positions <- sort(unique(exclude_positions))
  # Initialize the list to store subranges
  subranges <- list()
  # Add the first subrange from the beginning of the range to the first exclude position
  if (exclude_positions[1] > start(gr)) {
    subranges[[1]] <- GRanges(
      seqnames = seqnames(gr),
      ranges = IRanges(start = start(gr), end = exclude_positions[1] - 1)
    )
  }
  # Add subranges between the exclude positions
  for (i in seq_along(exclude_positions)[-1]) {
    subranges[[length(subranges) + 1]] <- GRanges(
      seqnames = seqnames(gr),
      ranges = IRanges(start = exclude_positions[i - 1] + 1, end = exclude_positions[i] - 1)
    )
  }
  # Add the last subrange from the last exclude position to the end of the range
  if (exclude_positions[length(exclude_positions)] < end(gr)) {
    subranges[[length(subranges) + 1]] <- GRanges(
      seqnames = seqnames(gr),
      ranges = IRanges(start = exclude_positions[length(exclude_positions)] + 1, end = end(gr))
    )
  }
  # Combine the subranges into a single GRanges object
  split_gr <- do.call(c, subranges)
  return(split_gr)
}

make_dist_matrix <- function(v, all_variants, distCutOff){
  
  n <- length(v)
  mat <- matrix(0, nrow = n, ncol = n)
  
  for (i in 1:n) {
    mat[i, ] <- abs(v - v[i])
  }
  rownames(mat) <- all_variants
  colnames(mat) <- all_variants
  mat[mat>distCutOff] <- 0
  mat[lower.tri(mat)] <- 0
  return(mat) 
}


get_general_info <- function(most_relevant_comb, all_comb, df_gene){
  func_start()
  if(str_count(paste(all_comb$status, collapse=' '),
               'diff')==nrow(all_comb)&
     nrow(all_comb)>=mean(as.numeric(df_gene$tcn))){ 
    
    info <- paste(
      'unclear: all affected is very likely because all',nrow(all_comb),
      'mutations are at different copies and tcn=', 
      mean(as.numeric(df_gene$tcn)),
      'tumor might be subclonal and therefore wt-copies left')
#  } else if(nrow(all_comb)>3){
    # info <-      paste('unclear: very high number of mutations...', 
    #                    'probably wt copies left but please check variant table')
  } else {
    # info <- paste0(most_relevant_comb$phasing, "-phasing of ", most_relevant_comb$comb, ": ", most_relevant_comb$status, 
    #                "; left wt cp: ", round(most_relevant_comb$wt_cp, 2)) 
    info <- NA
  }
  func_end()
  return(info)
}

get_unphased_info <- function(most_relevant_comb, all_comb){
  func_start()
  all_comb_unphased <- all_comb[which(all_comb$nstatus==0),]
  n_unphased <- nrow(all_comb_unphased)
  if(n_unphased>0){
    unphased <- TRUE
    unphased_info_pre <- paste(n_unphased, "of", 
                        nrow(all_comb), "unphased combinations:")
    if(most_relevant_comb$score==2){
      ## unphased ones do not matter
      unphased_info <- paste(unphased_info_pre, "But one affects all copies")  
    } else {
      unphased_info <- paste(unphased_info_pre, "Combs", 
                        paste(all_comb_unphased$comb, collapse = ", "), 
                        "might affect all copies") 
    }
  } else {
    unphased_info <- "all combinations solved"
    unphased <- FALSE
  }
  func_end()
  return(tibble(info=unphased_info, unphased=unphased))
}


eval_phasing_new <- function(all_comb, df_gene, printLog, verbose){
  func_start()
  
  rare_case <- ""
  most_relevant_comb <- all_comb %>%
    arrange(desc(score), wt_cp, desc(conf)) %>%
    .[1,]  
  
  if(most_relevant_comb$score==2){
    ## take most relevant combination
    gene_clear <- TRUE
    wt_cp_clear <- TRUE
    case <- "allaff"
    conf <- most_relevant_comb$conf
    conf_case <- "mostrel"
    eval <- most_relevant_comb$phasing
  } else {
    if(nrow(all_comb)==1){
      ## take the only combination
      eval <- most_relevant_comb$phasing
      if(most_relevant_comb$score==0){
        gene_clear <- FALSE
        wt_cp_clear <- FALSE
        case <- "onecomb_unphased"
        conf <- 0
        conf_case <- "poss_allaff_unphased"
      } else {
        gene_clear <- TRUE
        if(min(most_relevant_comb$min_poss_wt_cp, na.rm=TRUE)<0.5){
          conf <- most_relevant_comb$conf
          conf_case <- "poss_allaff_but_same"
        } else {
          conf <- 1
          conf_case <- "impos_allaff"
        }
        if(is.na(most_relevant_comb$wt_cp)){
          wt_cp_clear <- FALSE
          case <- "onecomb_wtNA"
        } else {
          wt_cp_clear <- TRUE
          case <- "onecomb_wtclear"
        }
      }
    } else {
      rare_case <- eval_rare_case(all_comb)
      ## evaluate all combinatioons
      if(0 %in% all_comb$score){
        ## potentially reachable all copies affectd by unphased combinations
        gene_clear <- FALSE
        wt_cp_clear <- FALSE
        case <- "morecomb_unphased"
        conf <- 0
        conf_case <- "poss_allaff_unphased"
        eval <- NA
      } else {
        ## each combination could be solved
        
        gene_clear <- TRUE
        if(any(!is.na(all_comb$wt_cp))){
          if(sum(!is.na(all_comb$wt_cp))==nrow(all_comb)){
            wt_cp_clear <- TRUE
            case <- "morecomb_allsolved"
          } else {
            if(min(all_comb$wt_cp, na.rm=TRUE)<=min(all_comb$min_poss_wt_cp, na.rm=TRUE)){
              ## if potential min wt copies are higher than clear ones
              wt_cp_clear <- TRUE
              case <- "morecomb_wtclear"
            } else {
              wt_cp_clear <- FALSE
              case <- "morecomb_wtNA"
            }            
          }

        } else {
          ## no wt copies annoatted 
          wt_cp_clear <- FALSE
          case <- "morecomb_onlyexc"
        }          
        
        if(min(all_comb$min_poss_wt_cp, na.rm=TRUE)<0.5){
            ## if this is the case, a possible diff would lead to all copies affected but it was phased as same
            ## for confidence measurement, we here need to take the phasing confidence which indicates this
            conf <- all_comb %>%
              mutate(p_comb=case_when(
                min_poss_wt_cp<0.5 ~ conf,
                TRUE ~ 1
              )) %>%
              pull(p_comb) %>%
              prod()
            rare_case <- " confidence adjusted"
            conf_case <- "morecomb_poss_allaff_but_same"
            
            all_phasing_annotations <- all_comb[which(!is.na(all_comb$phasing)),]
            if(nrow(all_phasing_annotations)>0){
              merged <- paste(all_phasing_annotations$phasing, collapse = ", ")
              eval <- case_when(
                str_detect(merged, "imbalance") ~ "imbalance",
                str_detect(merged, "haploblock") ~ "haploblock",
                str_detect(merged, "indirect") ~ "indirect",
                TRUE ~ "direct"
              )
            } else {
              eval <- "aff_cp"
            }
        } else {
            conf <- 1
            conf_case <- "morecomb_imposs_allaff"
            eval <- "aff_cp"
        }
      } ## each combination was solved
    } ## more than one comb
  } ## max score = 1
  
  if(gene_clear){
    gene_result <- most_relevant_comb %>%
      mutate(info=paste("most relevant phasing combination:",comb)) %>%
      select(gene, score, phasing, wt_cp) 
  } else {
    gene_result <- most_relevant_comb %>%
      select(gene) %>%
      mutate(info="unsolvable phasing combination", phasing=NA, score=0, wt_cp=NA)
  }
  wt_cp_range <- paste(round(min(all_comb$min_poss_wt_cp, na.rm=TRUE),2), "-",
        round(max(all_comb$max_poss_wt_cp, na.rm=TRUE),2))
  if(wt_cp_clear){
    wt_cp_result <- gene_result 
  } else {

    wt_cp_result <- gene_result %>%
      select(-wt_cp) %>%
      mutate(wt_cp=NA)
  }
  unplaus <- all_comb[which(!is.na(all_comb$unplausible)),]
  if(nrow(unplaus)>0){
  #if(any(!is.na(all_comb$unplausible))){
    if(any(str_detect(unplaus$unplausible, "1"))){
      
      issues <- unplaus %>%
        filter(str_detect(unplausible, "1")) %>%
        pull(comb) %>%
        paste(collapse = ", ") %>%
        paste("unplausible phasing combs:",.)
    } else {
      issues <- NA
    }    
  } else {
    issues <- NA
  }
  annotated_result <- wt_cp_result %>%
    mutate(conf=conf, 
           gene_clear=gene_clear, wt_cp_clear=wt_cp_clear, 
           case=paste0(case, rare_case), n_mut=nrow(df_gene),
           conf_case=conf_case,
           wt_cp_range=wt_cp_range,
           warning=issues,
           eval_by=eval
           )
  func_end()
  return(annotated_result)
}
#' @keywords internal
#' desicion tree if one variant already affects all copies in pre evaluation
eval_one_mut_affects_all <- function(df_gene, printLog, verbose){
  func_start(verbose)
  concern_info <- df_gene %>% 
    filter(vn_status==2) %>%
    .[1,] %>%
    mutate(n_mut=nrow(df_gene),
           info=paste("all copies affected by variant:", mut_id),
           conf=1,
           eval_by=ifelse(str_detect(class, "homdel"), "homdel","aff_cp")) %>%
    select(gene, n_mut, score=vn_status, conf, info, eval_by, wt_cp) %>%
    mutate(warning=NA, wt_cp_range=NA, phasing=NA)
  append_loglist(concern_info$info)
  func_end(verbose)
  return(concern_info)
}
#' @keywords internal
#' if no variant affecta all copies
eval_one_mut <- function(df_gene, printLog, verbose){
  func_start(verbose)
  concern_info <- df_gene %>%
    mutate(n_mut=1,
           conf=1,
           eval_by="aff_cp") %>%
    select(gene, n_mut, score=vn_status, conf, 
           info=pre_info, 
           wt_cp,
           eval_by) %>%
    mutate(warning=NA, wt_cp_range=NA, phasing=NA)
  append_loglist(concern_info$info)
  func_end(verbose)
  return(concern_info)
}
remove_duplicated_variants <- function(df_gene_raw, verbose){
  vm(as.character(sys.call()[1]), verbose, 1)
  n_vars_per_pos <- df_gene_raw %>%
    group_by(chr, pos) %>%
    tally() %>%
    filter(n>1)
  if(nrow(n_vars_per_pos)!=0){
    var_to_remove <- apply(n_vars_per_pos, 1, function(VAR){
      selected <- df_gene_raw %>% 
        filter(chr==VAR[["chr"]],
               pos==VAR[["pos"]]) %>%
        #filter(wt_cp==min(.$wt_cp)) %>%
        ## if both vars have the same impac we just select the first
        .[1,]
      warning(
        "more than one variant detected at position:", 
        VAR[["chr"]], ":", VAR[["pos"]],
        "\n  maybe an indel with shifted annotation of position?",
        "\nSelecting variant:",
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
  func_end(verbose)
  return(df_gene)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_replace_all str_count str_detect
#' @importFrom dplyr as_tibble group_by relocate bind_rows left_join select tally mutate filter nth relocate
predict_zygosity_genewise <- function(GENE, evaluation_per_variant, bamDna, 
                                      bamRna, showReadDetail,
                                      printLog, purity, sex, haploBlocks,
                                      phasedVcf, distCutOff, 
                                      verbose, logDir, refGen, somCna, 
                                      snpQualityCutOff, 
                                      phasingMode){
  #func_start(verbose)
  vm("predict_zygosity_genewise", verbose, 1) ## do not replcae!!!
  start_gene_eval <- Sys.time()
  increment_loglist()
  append_loglist(GENE)
  gene <- pre_info <-  chr <-  pos <- wt_cp <-  mut_id <- status <- . <- 
    score <- comb <- dist <- tcn <- info <- n <- all_comb <- phasing_info <- mat_phased_gene <- 
  mat_info_gene <- NULL 
  #start_gene_eval <- Sys.time()
  pre_df_gene <- evaluation_per_variant %>% filter(gene==GENE)    
  
  ## all germline variants which are lost in the tumor are excluded, 
  ## vn_status = -1
  df_gene_raw <- pre_df_gene%>%
    filter(vn_status>=0)
  ## check if variants at the same position are present
  df_gene <- remove_duplicated_variants(df_gene_raw, verbose)
  ## check if variants remain in that gene
  if(!nrow(df_gene)==0){
    if(any(df_gene$vn_status==2)){
      append_loglist("one of", nrow(df_gene), "variants affects all copies")
      eval_for_gene <- eval_one_mut_affects_all(df_gene, printLog, 
                                                verbose)
    ## (2): if only one variant is present  
    } else if(nrow(df_gene)==1){
      append_loglist("one variant that does ot affect all copies")
      eval_for_gene <- eval_one_mut(df_gene, printLog, 
                                    verbose)
    ## (3): more than one variant present in gene   
    } else {  
      append_loglist(nrow(df_gene),
            "heterozygous variants detected: Initializing haplotype phasing")
      #print(haploBlocks)
      ## returns two tibbles: first one direct phasing combinations
      ## second one indirect phasing combinations
      full_phasing_result <- phase(df_gene, bamDna, bamRna, showReadDetail,
                        purity, sex, haploBlocks, phasedVcf,
                        distCutOff, printLog, verbose, logDir, refGen, somCna, 
                        snpQualityCutOff, 
                        phasingMode)
      all_comb <- full_phasing_result[[1]] %>%
        mutate(gene=GENE)
      eval_for_gene <- eval_phasing_new(all_comb, df_gene,  
                                    printLog, verbose)
      phasing_info <- full_phasing_result[[2]] %>%
        mutate(gene=GENE)
      
      mat_phased_gene <- list()
      mat_info_gene <- list()
      mat_phased_gene[[GENE]] <- full_phasing_result[[3]]
      mat_info_gene[[GENE]] <- full_phasing_result[[4]]
    } 
    #print(eval_for_gene)
    end_gene_eval <- Sys.time()
    df_reduced <- eval_for_gene %>% as_tibble() %>%
      mutate(
        status=case_when(
          score==2 ~ "all_copies_affected",
          score==1 ~ "wt_copies_left",
          TRUE ~ "undefined"
        ),  
        eval_time_s=pt(start_gene_eval, end_gene_eval)
        ) 
    zygosity_gene <- list(pre_df_gene, 
                          list(df_reduced, 
                               all_comb,
                               phasing_info,
                               mat_phased_gene,
                               mat_info_gene)
                          )
  } else {
    ## all variants lost in tumor
    append_loglist("all variants lost in tumor")
    zygosity_gene <- list(pre_df_gene, NULL)
  } 
  log_message <- unlist(loglist) %>% paste(collapse = "\n")
  func_end(verbose)
  return(append(zygosity_gene, log_message))
}
#' @keywords internal
pt <- function(start_gene_eval, end_gene_eval){
  round(
    as.numeric(
      difftime(
        end_gene_eval, 
        start_gene_eval,
        units="secs")
    ),3)
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
                                             templateGenes, verbose){
  func_start(verbose)
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
  func_end(verbose)
  return(combined_uncovered)
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
set_global_variables <- function(debug, verbose, printLog){
  call_depth <<- 0
  log_depth <<- 0
  timelist  <<- list()
  if(debug==TRUE){
    dg <<- TRUE
  } else {
    dg <<- FALSE
  }
  if(verbose==TRUE){
    vb <<- TRUE
  } else {
    vb <<- FALSE
  }
  if(printLog==TRUE){
    plg <<- TRUE
  } else {
    plg <<- FALSE
  }
}
#' @keywords internal
catt <- function(printLog=FALSE, level, text){
  
  if(printLog==TRUE){
    message(rep("  ", level), text)
  }
}
increment_loglist <- function(){
  loglist <<- list()
  timelog <<- list()
}
append_loglist <- function(...){
  filler <- paste(rep("  ", call_depth), collapse = "")
  appendix <- paste(..., collapse = " ") %>%
    paste0(filler,.) %>%
    str_replace_all("\n", paste0("\n",filler))
  loglist <<- append(loglist, appendix)
  timelog <<- append(timelog, as.character(Sys.time()))
  if(plg==TRUE){
    message(appendix)
  }
}
add_timestamp <- function(){
  timelist[[call_depth]] <<- Sys.time()
}
remove_timestamp <- function(){
  timelist <<- timelist[c(1:call_depth)]
}
vm <- function(mes, verbose=FALSE, depth=0){
  if(depth>0){
    to_add <- "|-"
  }  
  if(depth<0){
    to_add <- "+-"
  }
  if(verbose==TRUE){
    to_print <- paste0(paste(rep("| ", call_depth),collapse=""), to_add,mes)
  } 
  if(depth>0){
    increment_call_depth()
    add_timestamp()
    timestamp <- ""
  }  
  if(depth<0){
    timestamp <- pt(timelist[[call_depth]], Sys.time()) %>%
      paste("  ~",.)
    decrement_call_depth()
    remove_timestamp()
  }
  if(verbose==TRUE){
    message(to_print, timestamp)
  }
}
