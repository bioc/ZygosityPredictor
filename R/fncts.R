calc_otsu <- function(data_vector){
  sorted_vector <- sort(data_vector)
  # Initialize variables to store best threshold and maximum variance
  best_threshold <- 0
  max_variance <- 0
  # Loop through all possible threshold indices (between sorted adjacent values)
  for (i in 1:(length(sorted_vector) - 1)) {
    #for (i in length(thds)-1) {
    # Calculate potential threshold as the middle value between adjacent values
    potential_threshold <- (sorted_vector[i]+sorted_vector[i+1])/2
    # Count the number of elements in each class
    num_lower <- i
    num_upper <- length(sorted_vector) - i
    # Calculate class probabilities
    w0 <- num_lower / length(sorted_vector)
    w1 <- num_upper / length(sorted_vector)
    # Calculate class means
    m0 <- sum(sorted_vector[1:i]) / num_lower
    m1 <- sum(sorted_vector[(i + 1):length(sorted_vector)]) / num_upper
    # Calculate intra-class variance
    intra_var <- w0 * w1 * (m0 - m1)^2
    # Update max_variance and best_threshold if necessary
    if (intra_var > max_variance) {
      max_variance <- intra_var
      best_threshold <- potential_threshold
    }
  }
  return(best_threshold)
}

prind <- function(x){
  if(dg==TRUE){
    #prind(1x)
  }
}
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
# prind <- function(x, debug=F){
#   if(debug==T){
#     #prind(1x)
#   }
# }
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

calc_left_wt_copies <- function(mtcn, nstatus, aff_copies1, aff_copies2){
  left_wt_copies <- ifelse(nstatus==2,
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
aggregate_probs <- function(ps){
  if(length(ps)>1){
    np <- ps[1]
    for (n in seq(1, length(ps)-1, 1)){
      np <- np+(1-np)*ps[n+1]
    }
    return(np)    
  } else {
    return(ps) 
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr as_tibble group_by mutate tibble tally 
#' @importFrom purrr set_names
classify_combination <- function(classified_reads, purity, printLog, 
                                 verbose=FALSE){
  func_start()
  result <- . <- fac <- NULL
  all_possible_results <- c('both', 'mut1', 'mut2', 'none', 
                            "dev_var", 'read_in_read', 
                            'skipped')
  ## calculate expected counts according to aff copies for both cases:
  ## same and diff
  exp_diff <- tibble(fac=all_possible_results[1:3],
                     exp=c(0,1,1)+0.0000001)
  exp_same <- tibble(fac=all_possible_results[1:3],
                     exp=c(1,0,0)+0.0000001)
  number <- classified_reads %>%
    mutate(fac=factor(result, levels = all_possible_results)) %>%
    group_by(fac, .drop = FALSE) %>%
    tally() %>% 
    column_to_rownames(var='fac') %>%
    t() %>%
    as_tibble()
#  both <- number[['both']]
 # mut1 <- number[['mut1']]
  #mut2 <- number[['mut2']]
  none_raw <- number[['none']]
  append_loglist(print_tibble(number))
  # append_loglist("both:", both, 
  #                "\nmut1:", mut1, 
  #                "\nmut2:", mut2, 
  #                "\nnone:", none_raw)
  ## calculate confidence from basecalls and mapping quality of read
  ## and aggregate them per classification result
  cr_conf <- classified_reads %>%
    rowwise() %>%
    mutate(pb1=ascii_to_dec(baseq1),
           pb2=ascii_to_dec(baseq2)) %>%
    ungroup() %>%
    mutate(
           mq1=10^(as.numeric(mapq1)/(-10)),
           mq2=10^(as.numeric(mapq2)/(-10)),
           ppos1=(1-pb1)*(1-mq1),
           ppos2=(1-pb2)*(1-mq2),
           p_result=sqrt(ppos1*ppos2)
           ) %>%
    mutate(fac=factor(result, levels = all_possible_results)) %>%
    group_by(fac, .drop = FALSE) %>%
    summarize(n=length(fac),
              prob_sum=sum(p_result))
  
  relevant_for_decision <- cr_conf %>%
    filter(fac %in% c("both","mut1", "mut2")) %>%
    select(fac, prob_sum)
  
  both <- relevant_for_decision %>%
    filter(fac=="both") %>%
    pull(prob_sum)
  mut1 <- relevant_for_decision %>%
    filter(fac=="mut1") %>%
    pull(prob_sum)
  mut2 <- relevant_for_decision %>%
    filter(fac=="mut2") %>%
    pull(prob_sum)
  

  
  append_loglist(print_tibble(relevant_for_decision))
  ## get total number of relevant classifications
  sum_rel <- sum(relevant_for_decision$prob_sum)
  if(sum_rel==0){
    status <- "null"
    p <- 1
    difference_norm_x <- ss_x <- sd_x <- ratio_norm_x <- NA
    nstatus <- evidence <- certainty <- confidence <- conf_log <- 0
    append_loglist("no evidence for any classification")
  } else {
    ## check for similarity of numbers in both, mut1 and mut2 by chi-squared
    sim_diff <- chisq.test(t(left_join(exp_diff, relevant_for_decision,
                                       by="fac") %>% 
                               select(exp, prob_sum)))
    sim_same <- chisq.test(t(left_join(exp_same, relevant_for_decision,
                                       by="fac") %>% 
                               select(exp, prob_sum)))
    ## extract x-squared and p value from chiquared
    ss_x=sim_same$statistic[[1]]
    sd_x=sim_diff$statistic[[1]]
    ## normalize x_squared by number of relevant classifications
    norm_ss_x <- ss_x/sum_rel
    norm_sd_x <- sd_x/sum_rel
    difference_norm_x <- abs(norm_ss_x-norm_sd_x)
    ratio_norm_x <- norm_sd_x/norm_ss_x    
    evid_diff <- case_when(
        #mut1+mut2==0 ~ 0,
        mut1==0|mut2==0 ~ 0.01,
        TRUE ~ mut1+mut2
    )
    evid_same <- case_when(
      both==0 ~ 0.01,
      TRUE ~2*both
    )
    if(sd_x<ss_x){
      status <- "diff"
      nstatus <- 2
      p <- sim_same$p.value
      evid_dec <- evid_diff
      evid_alt <- evid_same
    } else {
      status <- "same"
      nstatus <- 1
      p <- sim_diff$p.value
      evid_dec <- evid_same
      evid_alt <- evid_diff
    }  
    append_loglist("chisq against expected-diff_result: x_sqrd:", 
                   sd_x, "; p", sim_diff$p.value)
    append_loglist("chisq against expected-same_result: x_sqrd:", 
                   ss_x, "; p", sim_same$p.value)
    evidence <- case_when(
      evid_dec >= 20 ~ 3, #"high",
      evid_dec >= 3  ~ 2, #"medium"
      evid_dec > 1  ~ 1, #"low" 
      nstatus==1 ~ 1,
      TRUE ~ 0#"verylow"
    )
    cert_estimator <- evid_dec/evid_alt
    certainty <- case_when(
      cert_estimator >= 10 ~ 3,
      cert_estimator >= 5 ~ 2,
      cert_estimator >= 1 ~ 1,
      TRUE ~ 0
    )  
    conf_log <- log10((evid_dec^2)/(evid_alt*p))
    confidence <- case_when(
      conf_log<1 ~ 1,
      conf_log<3 ~ 2,
      conf_log<5 ~ 3,
      conf_log<10 ~ 4,
      TRUE ~ 5
    )
  }
  append_loglist("final status:", status)
  append_loglist("evidence:", evidence)
  append_loglist("certainty:", certainty)
  append_loglist("confidence:", confidence)
  status_table <- tibble(
    both=both,
    mut1=mut1,  
    mut2=mut2,
    dev_var=number[['dev_var']],
    no_overlap=sum(as.numeric(number[['read_in_read']]), 
                   as.numeric(number[['skipped']])),
    status=status,
    nstatus=nstatus,
    ss_x=ss_x,
    sd_x=sd_x,
    p=p,
    evid=evidence,
    cert=certainty,
    conf_log=conf_log,
    conf=confidence,
    diffnx=difference_norm_x,
    ratnx=ratio_norm_x,
    none_raw=none_raw,
    DNA_rds=nrow(classified_reads %>% filter(origin=="DNA")),
    RNA_rds=nrow(classified_reads %>% filter(origin=="RNA"))
    )
  func_end(verbose)
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
extract_subseq <- function(ref_pos, element_mut, parsed_read, 
                           length_indel, string){
  pos_in_seq <- ref_pos-element_mut$map_start+1
  subseq_in_element <- str_sub(element_mut[[string]], start=pos_in_seq, 
                               end=-1)
  subseq_in_following_elements <- 
    parsed_read[
      which(as.numeric(parsed_read$id)>as.numeric(element_mut$id)),] %>%
    pull(all_of(string)) %>%
    na.omit() %>%
    paste(collapse="")
  full_remaining <- paste0(subseq_in_element, subseq_in_following_elements)
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
        base <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
                               "seq")
        qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
                               "qual")
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
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "qual")
  } else {
    ## no inserttion detected.. extract base at position to get base quality
    info <- "no insertion detected: no I in cigar"
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "seq")
    qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "qual")
  }
  return(tibble(base=base, qual=qual, info=info, len_indel=len,
                exp_indel=expected_indel_detected))
}

cigar_element <- function(parsed_read, ref_pos){
  func_start()
  cigel <- parsed_read[which(between(rep(ref_pos, nrow(parsed_read)), 
                                     parsed_read$map_start, 
                                     parsed_read$map_end)),] %>%
    .[1,]
  func_end()
  return(cigel)
}

extract_base_at_refpos <- function(parsed_read, ref_pos, class, ref_alt, 
                                   ref_ref){
  func_start()
  ## extract elemnt containing reference position from parsed read
  element_mut <- cigar_element(parsed_read, ref_pos)
  ## if position is inside N or D type, the read does not cover the position
  ## N means for example skipped exon in RNA
  ## D wouldmean the position is deleted
  ## real deletions are still mapped in the M part
  if(element_mut$type %in% c("N", "D")){
    base_info <- tibble(base=NA, qual=NA, info="position skipped/deleted", 
                        len_indel=NA, exp_indel=NA)
  } else if(class %in% c("snv", "snp")){
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
  func_end()
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
    if(ref_class %in% c("snv", "snp")){
      pos_in_seq <- ref_pos-element_mut$map_start+1
      base <- 
        str_sub(element_mut$seq,
                start=pos_in_seq,
                end=pos_in_seq) 
      qual <- 
        str_sub(element_mut$qual,
                start=pos_in_seq,
                end=pos_in_seq)
      mut_at_read <- case_when(
        base==ref_alt ~ 1,
        base==ref_ref ~ 0,
        TRUE ~ -2
      ) 
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
  func_start()
  paired_reads <- bam[which(bam$qname==qname)] 
  read_start <- GenomicAlignments::start(paired_reads)
  seq <- as.character(paired_reads$seq)
  cigar <- paired_reads$cigar
  qual <- as.character(paired_reads$qual)
  mate <- c(1,2)
  ## parse cigar string according to query
  cigq <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar,
                                                        with.ops = T) 
  ## parse cigar string according to reference
  cigr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar,
                                                            with.ops = F) 
  raw_cigs <- lapply(mate, function(i){
    cigq_start <- GenomicAlignments::start(cigq[[i]])
    cigq_end <- GenomicAlignments::end(cigq[[i]])
    cigq_len <- length(cigq_start)
    cigr_start <- GenomicAlignments::start(cigr[[i]])
    cigr_end <- GenomicAlignments::end(cigr[[i]])
    seqq <- str_sub(seq[i],
                    start=cigq_start,
                    end=cigq_end) %>% na_if("")
    qualq <- str_sub(qual[i],
                     start=cigq_start,
                     end=cigq_end) %>% na_if("")
    width=pmax(GenomicAlignments::width(cigq[[i]]), 
               GenomicAlignments::width(cigr[[i]]))
    map_start=read_start[i]+cigr_end-width
    map_end=read_start[i]+cigr_end-1
    raw_cigs_new <- data.frame(
      width=width,
      type=names(cigq[[i]]),
      seq=seqq,
      qual=qualq,
      map_start=map_start,
      map_end=map_end,
      start=cigr_start,
      end=cigr_end,
      mate=mate[i],
      mapq=paired_reads$mapq[i],
      origin=paired_reads$origin[i]
    )
  }) %>%
    bind_rows() %>%
    mutate(id=c(1:nrow(.)))
  func_end()
  return(raw_cigs)
}
evaluate_base <- function(base_info, ref_alt, ref_ref){
  func_start()
  if(base_info$class %in% c("snv", "snp")){
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
  func_end()
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
  vm("core_tool", verbose, 1)
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
      mut1_in_read==1 ~ "mut1",
      mut2_in_read==1 ~ "mut2",
      TRUE ~ "dev_var"
    )    
  }
  func_end()
  return(c(qname=qname, result=final_assignment,
           origin=unique(parsed_read$origin),
           baseq1=base_info1$qual,
           mapq1=base_info1$mapq,
           baseq2=base_info2$qual,
           mapq2=base_info2$mapq
           )
         )
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
is_in_list <- function(row, my_list) {
  any(sapply(my_list, function(entry) all(row %in% entry)))
}
# recombine_three_combinations <- function(states){
#   #inp <- states
#   if(sum(length(which(states==0)), length(is.na(states)))==1&length(which(states==2))<2){
#     states[which(states==0)] <- max(states)
#   } 
#   #message(inp, " -> ", states)
#   return(states)  
# }
pick_next <- function(all_combinations_to_be_phased){
  all_combinations_to_be_phased %>%
    ungroup() %>%
    filter(prio==max(prio)) %>%
    .[1,] %>%
    pull(comb) %>%
    return()
}

recombine_three_combinations <- function(states){
  if(length(states)!=3){
    warning("THREE COMB SEMMS TO BE LONGER THAN 3!!!!!!!!!!")
  }
  inp <- states
  is_undefined <- which((is.na(states)|states==0))
  if(length(is_undefined)==1&length(which(states==2))<2){
    states[is_undefined] <- max(states, na.rm = TRUE)
  } 
  #message(inp, " -> ", states)
  return(states)  
}


subdivide_muts_in_three_combs <- function(all_comb){
  func_start()
  if(length(all_comb)==1){
    ret <- NULL
  } else {
    ret <- expand.grid(all_comb, all_comb, all_comb) %>% as_tibble() %>%
      filter(!(Var1==Var2|Var1==Var3|Var2==Var3)) %>%
      mutate_all(.funs = as.character) %>%
      rowwise() %>%
      mutate(ccid=paste(sort(c(Var1, Var2, Var3)), collapse = "_")) %>%
      select(ccid) %>%
      unique() %>%
      mutate(
        v1=unlist(str_split(ccid, "-|_")) %>% unique() %>% .[1],
        v2=unlist(str_split(ccid, "-|_")) %>% unique() %>% .[2],
        v3=unlist(str_split(ccid, "-|_")) %>% unique() %>% .[3],
        nv1=str_count(ccid, v1),
        nv2=str_count(ccid, v2),
        nv3=str_count(ccid, v3)) %>%
      filter(nv1==2, nv2==2, nv3==2) %>%
      pull(ccid)
  }
  func_end()
  return(ret)
}

define_next_priority <- function(master_combs, next_phasing){
  func_start()
  np <- master_combs %>%
    
    lapply(function(COMB){
     #print(COMB)
      tib <- COMB %>%
        str_split("_") %>%
        unlist() %>%
        as_tibble() %>%
        left_join(next_phasing, by=c("value"="comb")) 
      
      nna <- sum(is.na(tib$nstatus))
      nnull <- length(which(tib$nstatus==0))
      if(nnull>1){
        prio <- 0
      } else if(nna==1){
        prio <- 2
      } else {
        prio <- 1
      }
      
      if(any(str_count(tib$value, "m")>1)){
        prio <- prio+0.5
      }   
      
      tibf <- tib %>%
        mutate(prio=ifelse(is.na(nstatus), prio, 0),
               master_comb=COMB)
      
      return(select(tibf, comb=value, prio, master_comb))
      
    }) %>%
    bind_rows()  %>%
    left_join(next_phasing, by="comb")
  func_end()
  return(np)
}


#' @keywords internal
#' description follows
#' @importFrom stringr str_split
#' @importFrom purrr map_chr set_names
make_phasing_combinations <- function(df_mut_to_combine){
  func_start()
  . <- pos1 <- pos2 <- mut_id1 <- mut_id2 <- NULL
  comb <- outer(df_mut_to_combine$pos, df_mut_to_combine$pos, `paste`) %>%
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
  func_end()
  return(comb)
}
#' @keywords internal
#' description follows
check_for_overlapping_reads <- function(bamDna, bamRna,
                                        ref_chr1, 
                                        ref_chr2, 
                                        ref_pos1, 
                                        ref_pos2, 
                                        verbose=FALSE){
  #vm("checking DNA", verbose, 1)
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
    #vm("checking RNA", verbose, 1)
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

#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr tibble 
classify_reads <- function(mat_gene_relcomb, bamDna, bamRna, verbose){
  vm(as.character(sys.call()[1]), verbose, 1)
  #dist <- abs(as.numeric(line[['pos1']])-as.numeric(line[['pos2']]))
  
  ref_pos1 <- as.numeric(mat_gene_relcomb[,"pos"][[2]])  
  ref_pos2 <- as.numeric(mat_gene_relcomb[,"pos"][[1]])
  ref_chr1 <- as.numeric(mat_gene_relcomb[,"chr"][[2]])  
  ref_chr2 <- as.numeric(mat_gene_relcomb[,"chr"][[1]])
  
  ref_alt1 <- as.character(mat_gene_relcomb[,"alt"][[2]])
  ref_alt2 <- as.character(mat_gene_relcomb[,"alt"][[1]]) 
  
  ref_ref1 <- as.character(mat_gene_relcomb[,"ref"][[2]])
  ref_ref2 <- as.character(mat_gene_relcomb[,"ref"][[1]])
  
  ref_class1 <- as.character(mat_gene_relcomb[,"class"][[2]])
  ref_class2 <- as.character(mat_gene_relcomb[,"class"][[1]])
  #vm(line, verbose)
  #print(ref_class1)
  #print(ref_class2)
  bam <- check_for_overlapping_reads(bamDna,
                                     bamRna,
                                     ref_chr1,
                                     ref_chr2,
                                     ref_pos1,
                                     ref_pos2, 
                                     verbose)
  ##prind(1bam)
  ##prind(1ref_ref1)
  ##prind(1ref_ref2)
  if(!length(bam)==0){  
    #vm("reads detected, applying core function", verbose, 1)
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
    
    #vm("reads classified", verbose, 1)
    
  } else {
    #vm("no reads detected, return empty df", verbose, 1)
    classified_reads <- tibble()
  } 
  func_end(verbose)
  return(classified_reads)
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
      nstatus=0,
      score=1,
      conf=NA,
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

fill_phasing_status <- function(all_combinations, status, phasing_type, 
                                verbose){
  func_start(verbose)
  comb_template <- all_combinations %>%
      select(comb=comb_id)
  #prind(11)
  if(!is.null(status)){
    filled <- comb_template %>%
      left_join(status, by="comb")  
    #prind(12)
  } else {
    filled <- comb_template %>%
      mutate(status="null",
             nstatus=0,
             conf=0,
             phasing=phasing_type)
    #prind(13)
  }
  #prind(1filled)
  func_end(verbose)
  return(filled)
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

combine_main_muts_and_snps <- function(df_gene_filtered, snps){
  func_start()
  sub_comb <- bind_rows(
    df_gene_filtered %>%
      select(chr, pos, mut_id, ref, alt, class) %>%
      mutate_all(.funs = as.character),
    snps %>%
      mutate(class="snv") %>%
      select(chr=seqnames, pos=start, mut_id=snp_id, 
             ref=REF, alt=ALT, class) %>%
      mutate_all(.funs = as.character)
  ) %>%
    make_phasing_combinations() %>%
    ## THIS IS VERY IMPORTANT; FOR SOME REASON THE TRANSITION FROM FACTOR TO CHARACTERE IN THE NEXT FUNCTION CAUSES A -1 FOR THE CHROMOSOMSE!!!!!!
    mutate(chr1=as.character(chr1), chr2=as.character(chr2))
  func_end()
  return(sub_comb)
}




loadVcf <- function(vcf_in, chrom, region_to_load, refGen, verbose, somCna,which="all",
                    colname_gt="GT", colname_af="AF", colname_dp4="DP4", 
                    dkfz=FALSE){
  func_start(verbose)
  ## first check which format input vcf has
  #print(region_to_load)
  gr_list <- lapply(vcf_in, function(VCF){
    #vm(class(VCF), verbose)
    if(is(VCF, "TabixFile")){
      if(chrom %in% Rsamtools::seqnamesTabix(VCF)){
        ##prind(1"TabixFile")
        loadedVcf <- 
          readVcf(VCF, refGen, 
                  param=region_to_load)
        combVcf <- rowRanges(loadedVcf) #%>%
          #merge_sCNAs(.,somCna, verbose)
        if(length(combVcf)>0){
        ##prind(1loadedVcf)
          combVcf$ALT <- unlist(lapply(combVcf$ALT, 
                                       function(x){as.character(x[[1]][1])}))
          
          gt <- VariantAnnotation::geno(loadedVcf)[[colname_gt]][,1] %>% 
            as.character()
          af <- VariantAnnotation::info(loadedVcf)[[colname_af]]
          dp4 <- VariantAnnotation::info(loadedVcf)[[colname_dp4]] 
          
          combVcf$gt <- gt
          combVcf$dp4 <- dp4
          combVcf$af <- af
          vcf_info <- "variants_detected"
        } else {
          combVcf <- NULL
          vcf_info <- "no snps in roi"
        }
      } else {
        combVcf <- NULL
        vcf_info <- "vcf does not cover chromosome"
      }
    } else {
      if(str_detect(VCF, paste0("chr", chrom, "[^0-9]"))){
        ##prind(1"noTabix_det")
        ##prind("file detected")
        loadedVcf <- VariantAnnotation::readVcf(VCF, refGen)
        ##prind(loadedVcf)
        gt <- VariantAnnotation::geno(loadedVcf)[[colname_gt]] %>% 
          as.character()
        af <- VariantAnnotation::info(loadedVcf)[[colname_af]] 
        dp4 <- VariantAnnotation::info(loadedVcf)[[colname_dp4]]
        rangesVcf <- rowRanges(loadedVcf) 
        rangesVcf$gt <- gt
        rangesVcf$af <- af
        rangesVcf$dp4 <- dp4
        combVcf <-  subsetByOverlaps(rangesVcf, region_to_load)#)%>%
          
          #mergeByOverlaps(somCna) %>% as_tibble()
        if(length(combVcf)>0){
          combVcf$ALT <- unlist(lapply(combVcf$ALT, 
                                       function(x){as.character(x[[1]][1])}))
          vcf_info <- "variants_detected"
        } else {
          combVcf <- NULL
          vcf_info <- "no snps in roi"
        }
      } else {
        combVcf <- NULL
        vcf_info <- "no vcf file covering roi"
      }
    }
    ##prind(combVcf)
    return(list(combVcf, vcf_info))
  }) 
  ##prind("outside of loop")
  ##prind(gr_list)
  final_vcf <- lapply(gr_list, nth, 1) %>%
    compact() %>%
    Reduce(function(x,y)c(x,y),.)
  if(which=="HZ"){
    filtered_vcf <- final_vcf[which(final_vcf$gt %in% c("1/0", "0/1", "1|0", "0|1"))]
  } else {
    filtered_vcf <- final_vcf
  }
  func_end(verbose)
  return(filtered_vcf)
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
eval_haploblock_phasing <- function(do_haploblock_phasing,
                                    haploblock_phasing){
  if(do_haploblock_phasing){
    if(is.null(haploblock_phasing)){
      try <- TRUE
    } else {
      try <- FALSE
    }
  } else {
    try <- FALSE
  }
  return(try)
}

make_distance_tbl <- function(muts_in_hap, phased_snps_in_haploblock_ids, distCutOff){
  func_start()
  df_dist <- apply(muts_in_hap, 1, function(mut){
    phased_snps_in_haploblock_ids %>%
      mutate(dist=abs(start-as.numeric(mut[["pos"]])),
             mut_id=mut[["mut_id"]]) %>%
      return()
  }) %>% bind_rows() %>%
    filter(dist<distCutOff) %>%
    arrange(dist) %>%
    mutate(comb_id=paste(mut_id, snp_id, sep="-"))
  func_end()
  return(df_dist)
}


decide_following_phasing <- function(in_list, phasedVcf, verbose,
                                     phasing_mode="full", phasing_type, 
                                     haploBlocks=NULL){
  func_start()
  proceed <- FALSE
  file_presence <- FALSE
  conf_cutoff <- 3
  ## check if relevant files are present
  
   
  if(!is.null(phasedVcf)){
    if(phasing_type=="haploblock"){
      
      if(!is.null(haploBlocks)){
        file_presence <- TRUE
      }       
    } else {
      ## either imbalance or indirect
      file_presence <- TRUE
    }
  }
  if(file_presence==TRUE){
    if(phasing_mode=="full"){
      proceed <- TRUE
    } else {
      combined_status <- bind_rows(lapply(compact(in_list), nth, 1))
     #print(combined_status)
      if(length(unique(combined_status$comb))>1){
        if(any(combined_status$nstatus==0)){
          ## some are null
          proceed <- TRUE
        } else {
          ## no null results
          if(max(combined_status$conf)<conf_cutoff){
            ## only low confidence results
            proceed <- TRUE
          }
        }
      } else {
        ## only two muts
        max_status <- combined_status %>% filter(nstatus==max(nstatus))
        if(unique(max_status$nstatus)>0){
          if(max(max_status$conf)<conf_cutoff){
            ## only low confidence results
            proceed <- TRUE
          }
        } else {
          ## only null results
          proceed <- TRUE
        }
      }     
    }
  }
  func_end()
  return(proceed)
}


select_status <- function(status_direct, status_haploblock){
  if("diff" %in% c(status_direct, status_haploblock)){
    return("diff")
  } else if("same" %in% c(status_direct, status_haploblock)){
    return("same")
  } else {
    return("null")
  }
}
print_tibble <- function(tbl_in){
  res <- paste(as.character(knitr::kable(tbl_in)), collapse="\n")
  return(res)
}
aggregate_phasing <- function(in_list, all_combinations, geneDir, verbose){
  func_start(verbose)
  phasing_types <- factor(c("direct","indirect", "haploblock", "imbalance"))
  corr_factors <- tibble(phasing=phasing_types,
                         factor=c(1,1,1,0.5))
  all_phasing_types <- lapply(in_list, nth, 1) %>% bind_rows() %>%
    left_join(corr_factors, by="phasing") %>%
    mutate(conf_type=conf*factor)
  append_loglist(print_tibble(all_phasing_types %>% 
                                mutate(gene=basename(geneDir))))
  store_log(geneDir, all_phasing_types %>% 
              mutate(gene=basename(geneDir)), "full_phasing_status.tsv")
  per_comb <- lapply(unique(all_phasing_types$comb), function(COMB){
    rel_comb <- all_phasing_types %>%
      filter(comb==COMB) 
    if(1 %in% rel_comb$nstatus&2 %in% rel_comb$nstatus){
      warning("different phasing approaches assigned different states")
      final_result <- rel_comb %>%
        arrange(desc(conf_type), phasing) %>%
        .[1,]
    } else {
      final_result <- rel_comb %>%
        arrange(desc(nstatus), desc(conf_type), phasing)%>%
        .[1,]
    }
    return(final_result)
  }) %>%
    bind_rows() %>%
    left_join(all_combinations, by=c("comb"="comb_id_sorted")) %>%
    rowwise() %>%
    mutate(
    tcn=min(tcn1, tcn2),
    wt_cp=calc_left_wt_copies(tcn,
                              nstatus,
                              aff_cp1,
                              aff_cp2),
    score=case_when(
      nstatus==2&wt_cp<0.5 ~ 2,
      nstatus==0 ~ 0,
      TRUE ~ 1),
    ) %>% 
    select(1:3, phasing, conf=conf_type, wt_cp, score)
  
  func_end(verbose)
  return(per_comb)
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
get_string_status <- function(nstatus){
  case_when(nstatus==2 ~ "diff",
            nstatus==1 ~ "same",
            TRUE ~ "null"
  )
}
load_snps <- function(df_gene, phasedVcf, refGen, distCutOff, verbose, somCna, which="all"){
  func_start()
  region_to_load <- paste0(unique(df_gene$chr), ":", min(df_gene$pos)-distCutOff,
                           "-", max(df_gene$pos)+distCutOff) %>%
    GRanges()  %>%
    split_genomic_range(.,df_gene$pos)
  
  loaded_vcf_hz <- loadVcf(phasedVcf, unique(df_gene$chr), region_to_load, refGen, 
                           verbose, somCna)#, which=which) 
  #define_region_to_load(df_gene)
  snps <- loaded_vcf_hz %>% as_tibble() %>%
    mutate(snp_id=paste0("s", c(1:nrow(.))+nrow(df_gene)))
  
  if(!"af" %in% names(snps)){
    if("dp4" %in% names(snps)){
      snps$af <- lapply(seq(1,nrow(snps)), function(i){
        sum(snps$dp4[[i]][c(3,4)])/sum(snps$dp4[[i]])
      }) %>% as.numeric()  
    } 
  } 
  
  
  gr_snps <- GRanges(snps %>% select(1:3, snp_id))
  merged_haploblocks <- mergeByOverlaps(gr_snps, haploBlocks) %>%
    as_tibble() %>%
    select(all_of(names(.) %>% .[which(!str_detect(.,"\\."))]))
  merged_somCna <- mergeByOverlaps(gr_snps, somCna) %>%
    as_tibble() %>%
    select(all_of(names(.) %>% .[which(!str_detect(.,"\\."))]))
  
  annotated_snps <- left_join(
    snps,
    merged_haploblocks,
    by="snp_id"
  ) %>%
    left_join(
      .,
      merged_somCna,
      by="snp_id"
    )
  
  func_end()
  return(annotated_snps)
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
append_matrices <- function(df_status, iterate=TRUE){
  mat_old <- mat_phased
  lapply(seq(1,nrow(df_status)), function(C){
    #split_list <- str_split(df_status[C,]$comb, "-") %>% unlist()
    #print(split_list)
    #mat_phased[split_list[1], split_list[2]] <<- df_status[C,]$nstatus
    #mat_conf[split_list[1], split_list[2]] <<- df_status[C,]$conf
    mat_phased[df_status[C,]$comb] <<- df_status[C,]$nstatus
    mat_conf[df_status[C,]$comb] <<- df_status[C,]$conf
    mat_info[df_status[C,]$comb] <<- df_status[C,]$phasing
  })
  if(iterate==TRUE){
    something_changed <- sum(mat_old, na.rm = T)!=sum(mat_phased, na.rm = T)
    mat_new <- mat_phased
    mat_new[is.na(mat_new)] <- 0
    mat_old[is.na(mat_old)] <- 0
    changes <- which(!mat_new==mat_old) 
    i <- 1
    ##solved_master_combs <<- list()
    while(i<=length(changes)){
      message("left changes:", changes)
      ## something changed
      ##print("something changed")
      ## get mut/snp ids that have changed
      m12 <- get_xy_index(changes[i], nrow(mat_new))
      status <- mat_new[changes[i]]
      ## get third mut/snp id that also has a non-null connection to one of the first ones 
      m1_same_conn <- c(m12[1], which(mat_new[m12[1],]==1),
                 which(mat_new[,m12[1]]==1)) %>% sort()
        
      m2_same_conn <- c(m12[2],
                 which(mat_new[m12[2],]==1),
                 which(mat_new[,m12[2]]==1)) %>% sort()
        
      if(sum(c(length(m1_same_conn), length(m2_same_conn)))>2){
        if(status==1){
          print("state 1")
          conns <- sort(unique(c(m1_same_conn, m2_same_conn)))
              for(i in 1:(length(conns)-1)){
                #print(i)
                for(j in (i+1):length(conns)){
                 # print(j)
                  #print(paste(conns[i],conns[j], sep="-"))
                  if(is.na(mat_phased[conns[i],conns[j]])|mat_phased[conns[i],conns[j]]==0){
                    mat_phased[conns[i],conns[j]] <<- 1
                    mat_info[conns[i],conns[j]] <<- changes[i]
                  }
                  
                }
              }
        } else {
          print("state 2")
          ## if status == diff
          for(i in m1_same_conn){
            #print(i)
            for(j in m2_same_conn){
              # print(j)
              #print(paste(conns[i],conns[j], sep="-"))
              sv <- sort(c(i,j))
              if(is.na(mat_phased[sv[1], sv[2]])|mat_phased[sv[1], sv[2]]==0){
                mat_phased[sv[1], sv[2]] <<- 2
                mat_info[sv[1], sv[2]] <<- changes[i]
              }
            }
          }
        }      
      }
  
      ## first annoatte all vars that are "same" in matrices

      # 
      # all_possible_changes <- sort(unique(c(m12, m3)))
      # 
      # 
      # 
      # pos_master_combs <- expand.grid(all_possible_changes, all_possible_changes, all_possible_changes) %>%
      #   .[which(!(.$Var1==.$Var2|.$Var1==.$Var3|.$Var2==.$Var3)),] %>%
      #   apply(.,1,sort) %>% t() %>% unique()
      # 
      # filtered_combs <- pos_master_combs[!apply(pos_master_combs, 1, is_in_list, solved_master_combs), ] 
      # 
      # 
      # 
      # if(length(filtered_combs)>0){
      #   if(is.matrix(filtered_combs)==FALSE){
      #     filtered_combs <- matrix(filtered_combs, ncol=3)
      #   }
      #   new_changes <- apply(filtered_combs, 1, function(I){
      #     rel_combs <- mat_phased[I,I]  #%>% recombine_three_combinations() %>%
      #     print(rel_combs)
      #     recomb_combs <- recombine_three_combinations(rel_combs%>% .[upper.tri(.)])
      #     
      #     rel_combs[upper.tri(rel_combs)] <- recomb_combs
      #     
      #     ## now if sth changed, put it in global matrix
      #     if(!0 %in% recomb_combs){
      #       solved_master_combs <<- append(solved_master_combs, list(I))
      #     }
      #     print(recomb_combs)
      #     print(rel_combs)
      #     print(mat_phased[I,I])
      #     print(rel_combs!=mat_phased[I,I])
      #     if(isTRUE(any(rel_combs!=mat_phased[I,I]))){
      #       n_changed <- which(rel_combs!=mat_phased[I,I])
      #       coords <- get_xy_index(n_changed, 3)
      #       changex <- I[coords[1]]
      #       changey <- I[coords[2]]
      #       new_val <- rel_combs[n_changed]
      #       new_conf <- mat_conf[I,I] %>% .[which(.!=0)] %>% mean()
      #       mat_phased[changey, changex] <<- new_val
      #       mat_conf[changey, changex] <<- new_conf
      #       mat_info[changey, changex] <<- paste(I, collapse="-")
      #       return(get_single_index(changex, changey, nrow(mat_phased)))
      #     } else {
      #       return(NULL)
      #     }
      #   }) %>% compact() %>%
      #     unlist()
      #   changes <- c(changes, new_changes) %>% unique()
      # }
      # 
      # if(length(changes)==1){
      #   something_changed <- FALSE
      # } else {
      #   changes <- changes[2:length(changes)]
      # }  
      i <- i+1 
    } 
  }
  
  #print(mat_phased)
}
get_single_index <- function(x,y,nr){
  
  nr*(x-1)+y
  
}

get_xy_index <- function(inp, nr){
  lapply(inp, function(i){
    y <- i%%nr
    if(y==0){
      y_ret <- nr
      x <- (i-y)/nr 
    } else {
      y_ret <- y
      x <- (i-y)/nr+1
    }
    return(c(x=x, y=y_ret))  
  }) %>% Reduce(function(x,y)rbind(x,y),.)
  
}
create_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  mat_dist <<- make_dist_matrix(all_pos, all_variants, distCutOff) 
  #if(sum(mat_dist[seq(1,length(df_gene$mut_id)),])>0){
  main_muts <<- all_variants
  main_pos <<- all_pos
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA    
  mat_conf <<- mat_phased
  mat_info <<- mat_phased
  #} else {
    ## no snps closer than distCutoff to main muts
  #}
}
append_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  
  mat_dist <<- make_dist_matrix(c(main_pos, all_pos), 
                               c(main_muts, all_variants), 
                               distCutOff) 
  mat_phased_main <- mat_phased
  mat_conf_main <- mat_conf
  mat_info_main <- mat_info
    
  
  
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA 
  mat_conf <<- mat_phased
  mat_info <<- mat_phased
  #print(mat_phased)
  
  nm <- length(main_muts)
  #print(mat_phased[nm, nm])
  #print(mat_phased_main)
  mat_phased[c(1:nm), c(1:nm)] <<- mat_phased_main
  mat_conf[c(1:nm), c(1:nm)] <<- mat_conf_main
  mat_info[c(1:nm), c(1:nm)] <<- mat_info_main
  
}
prioritize_combination <- function(mut_ids, snp_ids){
  unrel_rows <- rownames(mat_phased) %>% .[which(!. %in% c(mut_ids, snp_ids))]
  mat_tmp <- mat_dist
  mat_tmp[unrel_rows,] <- NA
  mat_tmp[,unrel_rows] <- NA
  mat_tmp[c((length(main_muts)+1):nrow(mat_phased)),] <- NA
  unphased <- which(is.na(mat_phased))
  relevant_dists <- mat_tmp[unphased] 
  relevant_dists[which(relevant_dists==0)] <- NA
  shortest <- unphased[which(relevant_dists==min(relevant_dists, na.rm=T))]
  next_comb <- get_xy_index(shortest, nrow(mat_phased))
  return(next_comb)
}
unphased_main <- function(){
  lm <- seq(1,length(main_muts),1)
  which(is.na(mat_phased[lm, lm]))
}
unknown_main <- function(){
  lm <- seq(1,length(main_muts),1)
  mat_tmp <- mat_phased[lm, lm]
  which(upper.tri(mat_tmp)&mat_tmp==0)
}
add_snps_to_matrices <- function(snps){
  snps_hap <- snps[which(!is.na(snps$block_final)),] %>%
    select(mut_id, block_final, gt=gt_final)
  known_combs <- lapply(unique(snps_hap$block_final) %>% .[which(.!=0)], function(HAP_ID){
    
    snps_in_hap <- snps_hap %>%
      filter(block_final==HAP_ID) %>%
      mutate(id=as.numeric(str_match(mut_id, "\\d+")))
    res_list <- list()
    
    for (i in 1:(nrow(snps_in_hap)-1)) {
      result_matrix <- matrix(0, nrow = nrow(snps_in_hap)-i, ncol = 2)
      for (j in (i+1):nrow(snps_in_hap)) {
        comparison_result <- ifelse(snps_in_hap$gt[i] == snps_in_hap$gt[j], 1, 2)
        # Store the results in the result_matrix
        result_matrix[nrow(snps_in_hap)-j+1, 1] <- get_single_index(snps_in_hap$id[j], snps_in_hap$id[i], nrow(mat_phased))
        result_matrix[nrow(snps_in_hap)-j+1, 2] <- comparison_result
        
      }
      res_list[[i]] <- result_matrix
    }
    
    to_add <- Reduce(function(x,y)rbind(x,y),res_list) #%>% t() %>%
    colnames(to_add) <- c("comb", "nstatus") 
    new <- to_add %>%
      as_tibble() %>%
      mutate(conf=5, phasing=0, hap_id=HAP_ID) #%>
    
    return(new)
  })  %>% bind_rows()
  
  append_matrices(known_combs, iterate=FALSE)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
phase <- function(df_gene, bamDna, bamRna, 
                  showReadDetail, purity, sex, haploBlocks, phasedVcf,
                  distCutOff, printLog, verbose, logDir, refGen, somCna, 
                  snpQualityCutOff, 
                  phasingMode){
  vm(as.character(sys.call()[1]), verbose, 1)
  #haploblock_phasing <- NULL
  if(!is.null(logDir)){
    geneDir <- file.path(logDir, unique(df_gene$gene))
    dir.create(geneDir)
  }
  ##prind(1df_gene)
  ## (1): define all combinations of variants to be phased
  #all_combinations <- make_phasing_combinations(df_gene) 
  create_phasing_matrices(df_gene$mut_id, df_gene$pos, distCutOff)
  solved_master_combs <<- list()
  mat_gene <- as.matrix(df_gene)
  rownames(mat_gene) <- mat_gene[,"mut_id"]
  
  unphased <- unphased_main()
  
  
  append_loglist(length(unphased), "main combinations to phase,", 
                 abs(length(unphased)-length(mat_phased[upper.tri(mat_phased)])), 
                 "are over distCutOff")
  ## (2): perform direct phasing between variants (read-based)
  direct_phasing <- perform_direct_phasing(mat_gene, bamDna, bamRna, 
                                           purity, verbose, printLog, 
                                           showReadDetail, geneDir)
  
  if(length(unknown_main())>0&!is.null(phasedVcf)){
    ## missing combinations in main muts --> start secondary phasing approaches
    lsnps <- load_snps(df_gene, phasedVcf, refGen, distCutOff, verbose, somCna) %>%
      dplyr::rename(mut_id=snp_id)
    snps <- lsnps %>%
      rowwise() %>%
      mutate(
        aff_cp=aff_germ_copies(seqnames, af, tcn, purity, sex)
      ) %>%
      ungroup() %>%
      mutate(
        gt_cna_min=str_split(gt_cna, ":") %>% unlist() %>% min(),
        gt_cna_max=str_split(gt_cna, ":") %>% unlist() %>% max(),
        gt_seg=case_when(
          all_imb==FALSE ~ NA,
          round(aff_cp) == gt_cna_min ~ "0|1",
          round(aff_cp) == gt_cna_max ~ "1|0",
        ),
        gt_final=case_when(
          gt=="1|1"|gt=="0|0" ~ NA,
          !is.na(hap_id)&str_detect(gt, "\\|") ~ gt,
          !is.na(gt_seg) ~ gt_seg,
          TRUE ~ NA
        ),
        block_final=case_when(
          gt=="1|1"|gt=="0|0" ~ NA,
          !is.na(hap_id)&str_detect(gt, "\\|") ~ paste0("h", hap_id),
          !is.na(gt_seg) ~ paste0("s", seg_id),
          TRUE ~ NA
          
        )
      )
    
    append_phasing_matrices(snps$mut_id, snps$start, distCutOff)
    
    add_snps_to_matrices(snps)
    
    ## build main daatframe of muts and snps and annotate haploblocks and segments of allelic imbalance
    df <- bind_rows(
      df_gene %>% select(chr, pos, ref, alt, class, af, mut_id),
      snps %>% select(chr=seqnames, pos=start, ref=REF, alt=ALT, af, mut_id, gt=gt_final) %>% mutate(class="snp")
    ) 
    
    phasing_info <- tibble()
    to_phase <- prioritize_combination(main_muts, 
                                       snps$mut_id)          
    i <- 1
    while(!is.null(to_phase)){
      print(to_phase)
      
      mat_gene_relcomb <- df[to_phase,] %>% as.matrix() 
      comb <- get_single_index(to_phase[1], 
                               to_phase[2], 
                               nrow(mat_phased))
      
      classified_main_comb <- phase_combination(mat_gene_relcomb, comb, 
                                                bamDna, bamRna, verbose, geneDir, 
                                                geneDir, comb, showReadDetail)  
      phasing_info <- bind_rows(phasing_info,
                                classified_main_comb)
      to_phase <- prioritize_combination(main_muts, 
                                         snps$mut_id) 
      i <- i +1
    }
    print(mat_phased)
    stop("stop here")
    
   } ## secondary phasing done
  
  ## reconstruct phasing results
  
  
  
  
  


 
  
  # rownames(mat) <- all_variants
  # colnames(mat) <- all_variants
  # mat[mat>distCutOff] <- 0
  # mat[lower.tri(mat)] <- 0
  # region_to_load <- paste0(unique(df_gene$chr), ":", min(df_gene$pos)-distCutOff,
  #                          "-", max(df_gene$pos)+distCutOff) %>%
  #   GRanges()  %>%
  #   split_genomic_range(.,df_gene$pos)
  # 
  # loaded_vcf_hz <- loadVcf(phasedVcf, unique(df_gene$chr), region_to_load, refGen, 
  #                          verbose, "HZ") 
  # #define_region_to_load(df_gene)
  # snps <- loaded_vcf_hz %>% as_tibble() %>%
  #   mutate(snp_id=paste0("s", c(1:nrow(.))+nrow(df_gene)))
  
  
  # mat_main_comb <- 
  # 
  # full_snps_in_gene <- load_all_snps_in_relevant_gene_region(df_gene, phasedVcf, 
  #                                       refGen, verbose, 
  #                                       distCutOff)
  ## load relevant snps for haploblock and imbalance phasing
  # indirect_phasing <- perform_indirect_phasing(df_gene, direct_phasing, phasedVcf, 
  #                                              refGen, verbose, phasingMode, 
  #                                              all_combinations, showReadDetail, 
  #                                              geneDir, distCutOff)  
  

  imbalance_phasing <- perform_imbalance_phasing(all_combinations, direct_phasing, 
                                                 haploblock_phasing, df_gene, 
                                                 somCna, phasedVcf, bamDna, bamRna, 
                                                 purity, sex, distCutOff, 
                                                 printLog, 
                                                 verbose, geneDir, refGen, snpQualityCutOff, 
                                                 phasingMode)  
  

  

  phasing_results_list <- list(direct_phasing,
                               #indirect_phasing,
                               haploblock_phasing,
                               imbalance_phasing)
  aggregated_phasing <- aggregate_phasing(phasing_results_list, 
                                          all_combinations, geneDir, verbose)
 #print(aggregated_phasing)
  func_end()
  return(aggregated_phasing) 
}


eval_phasing_new <- function(all_comb, df_gene, printLog, verbose){
  func_start()
  req_output_cols <- c("gene", "n_mut", "score", "conf",  "info")
  if(nrow(all_comb)==1){
    ## only one combination
    gene_phasing <- all_comb %>%
      mutate(info=paste0(phasing, "-phasing of ", comb, ": ", status, 
                         "; left wt cp: ", round(wt_cp, 2)))
  } else {
    
    ## check if combinationa re missing and try to fill them
    
    if(max(all_comb$score)==2){  
      ## at least two combinations one of whcih affects all copies
      gene_phasing <- all_comb %>% 
        filter(wt_cp==min(all_comb$wt_cp)) %>%
        .[1,] %>%
        mutate(info=paste0(phasing, "-phasing of ", comb, ": ", status, 
                           "; left wt cp: ", round(wt_cp, 2))) 
      
      
    } else if(str_count(paste(all_comb$status, collapse=' '),
                        'diff')==nrow(all_comb)&
              nrow(all_comb)>=mean(as.numeric(df_gene$tcn))){ 
      
      gene_phasing <- all_comb %>% 
        filter(wt_cp==min(all_comb$wt_cp)) %>%
        .[1,]%>% 
        mutate(
          info=paste(
          'unclear: all affected is very likely because all',nrow(all_comb),
          'mutations are at different copies and tcn=', 
          mean(as.numeric(df_gene$tcn)),
          'tumor might be subclonal and therefore wt-copies left'))
      
      
      
    } else if(nrow(all_comb)>3){
      
      gene_phasing <- all_comb %>% 
        filter(wt_cp==min(all_comb$wt_cp)) %>%
        .[1,] %>%
        mutate(info=      paste('unclear: very high number of mutations...', 
              'probably wt copies left but please check variant table'))
      
    } else {
      gene_phasing <- all_comb %>% 
        filter(wt_cp==min(all_comb$wt_cp)) %>%
        .[1,] %>%
        mutate(info="test")
    }
    
    
    
  }
    
    
      
  func_end()
  return(gene_phasing %>%
           mutate(n_mut=nrow(df_gene)) %>%
           select(all_of(req_output_cols)))
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
           conf=NA) %>%
    select(gene, n_mut, score=vn_status, conf, info)
  append_loglist(concern_info$info)
    #filter(str_detect(pre_info, "all copies affected")) %>% .[1,]
  func_end(verbose)
  return(concern_info)
}
#' @keywords internal
#' if no variant affecta all copies
eval_one_mut <- function(df_gene, printLog, verbose){
  func_start(verbose)
  concern_info <- df_gene %>%
    mutate(n_mut=1,
           conf=NA) %>%
    select(gene, n_mut, score=vn_status, conf, info=pre_info)
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
    score <- comb <- dist <- tcn <- info <- n <- all_comb <- NULL 
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
      
      ## returns two tibbles: first one direct phasing combinations
      ## second one indirect phasing combinations
      full_phasing_result <- phase(df_gene, bamDna, bamRna, showReadDetail,
                        purity, sex, haploBlocks, phasedVcf,
                        distCutOff, printLog, verbose, logDir, refGen, somCna, 
                        snpQualityCutOff, 
                        phasingMode)
      all_comb <- full_phasing_result %>%
        mutate(gene=GENE)
      eval_for_gene <- eval_phasing_new(all_comb, df_gene,  
                                    printLog, verbose)
    } 
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
                               all_comb)
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
