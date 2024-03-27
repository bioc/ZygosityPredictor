func_end <- function(ZP_env){
  vm("end", -1, ZP_env)
}
func_start <- function(ZP_env){
  mes <- paste0(as.character(sys.call(-1)[1]))
  vm(mes, 1, ZP_env)
}
increment_log_depth <- function(ZP_env) {
  ZP_env$global_ZygosityPredictor_variable_debug_log_depth <- ZP_env$global_ZygosityPredictor_variable_debug_log_depth + 1
}
decrement_log_depth <- function(ZP_env) {
  ZP_env$global_ZygosityPredictor_variable_debug_log_depth <- ZP_env$global_ZygosityPredictor_variable_debug_log_depth - 1
}
increment_call_depth <- function(ZP_env) {
  ZP_env$global_ZygosityPredictor_variable_debug_call_depth <- ZP_env$global_ZygosityPredictor_variable_debug_call_depth + 1
}
decrement_call_depth <- function(ZP_env) {
  ZP_env$global_ZygosityPredictor_variable_debug_call_depth <- ZP_env$global_ZygosityPredictor_variable_debug_call_depth - 1
}
#' @importFrom readr write_tsv
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
#' @importFrom purrr compact
#' @importFrom dplyr summarize pull
insert_missing_cnv_regions <- function(somCna, geneModel, sex, ploidy, ZP_env){
  vm(as.character(sys.call()[1]), 1, ZP_env)
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
  func_end(ZP_env)
  return(comb_somCna)
}

#' @keywords internal
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate select mutate_at rowwise
#' @importFrom stringr str_detect
merge_sCNAs <- function(obj, somCna, ZP_env){
  func_start(ZP_env)
  alt <- af <- tcn <- tcn_assumed <- cna_type <- gene <- ref <- NULL
  merged <- obj %>%
    mergeByOverlaps(somCna) %>% as_tibble() %>%
    mutate(cna_type=ifelse(str_detect(cna_type, "LOH"), "LOH", "HZ")) %>%
    select(chr=1, pos=2, gene, ref, alt, af, tcn, cna_type, all_imb, gt_cna, 
           seg_id, tcn_assumed) %>%
    mutate_at(.vars = c("af", "tcn"), .funs=as.numeric) %>%
    rowwise() #%>%
  func_end(ZP_env)
  return(merged)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate filter
prepare_somatic_variant_table <- function(somSmallVars, templateGenes, 
                                          somCna, purity, sex, ZP_env){
  func_start(ZP_env)
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- 
    aff_cp <- wt_cp <- . <- tcn_assumed <- NULL
  if(!is.null(somSmallVars)){
    reduced_to_templateGenes <- 
      somSmallVars[which(somSmallVars$gene %in% templateGenes)]
    if(length(reduced_to_templateGenes)>0){
      tbl_prepared_variants <- reduced_to_templateGenes %>%
        merge_sCNAs(., somCna, ZP_env) %>%
        
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
  func_end(ZP_env)
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
                                       byTcn, sex, ploidy, include_dels, ZP_env){
  #print("jou")
  func_start(ZP_env)
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
                   gt_cna=NA,
                   seg_id=NA,
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
  func_end(ZP_env)
  return(full_CNAs)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom dplyr mutate select
prepare_germline_variants <- function(germSmallVars, somCna, purity, sex, ZP_env){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- aff_cp <- 
    origin <- pos <- wt_cp <- pre_info <- . <- tcn_assumed <- NULL
  func_start(ZP_env)
  if(is.null(germSmallVars)){
    df_germ <- NULL
  } else {
    df_germ <- germSmallVars %>% 
      merge_sCNAs(., somCna, ZP_env) %>%
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
             gt_cna,
             seg_id,
             aff_cp,
             wt_cp,
             pre_info,
             vn_status
      ) 
  }
  func_end(ZP_env)
  return(df_germ)
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

#' @importFrom knitr kable
print_tibble <- function(tbl_in){
  # res <- paste(as.character(knitr::kable(tbl_in)), collapse="\n")
  res <- paste(as.character(knitr::kable(tbl_in)), collapse="\n")
  return(res)
}
#' @importFrom GenomicRanges GRanges seqnames start end
#' @importFrom IRanges IRanges
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
#' @importFrom stringr str_detect %>%
#' @importFrom dplyr filter pull
find_unplaus_issues <- function(unplaus){
  if(nrow(unplaus)>0){
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
}
#' @importFrom dplyr arrange case_when desc mutate pull select
#' @importFrom stringr %>% str_detect 
eval_phasing_new <- function(all_comb, df_gene, printLog, ZP_env){
  func_start(ZP_env)
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
        gene_clear <- wt_cp_clear <- FALSE
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
                str_detect(merged, "imbalance") ~ "imbalance-phasing",
                str_detect(merged, "haploblock") ~ "haploblock-phasing",
                str_detect(merged, "indirect") ~ "indirect-phasing",
                TRUE ~ "direct-phasing"
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
  issues <- find_unplaus_issues(unplaus)
  annotated_result <- wt_cp_result %>%
    mutate(conf=conf, 
           gene_clear=gene_clear, wt_cp_clear=wt_cp_clear, 
           case=paste0(case, rare_case), n_mut=nrow(df_gene),
           conf_case=conf_case,
           wt_cp_range=wt_cp_range,
           warning=issues,
           eval_by=eval
           )
  func_end(ZP_env)
  return(annotated_result)
}
#' @keywords internal
#' desicion tree if one variant already affects all copies in pre evaluation
#' @importFrom dplyr mutate select filter
#' @importFrom stringr str_detect
eval_one_mut_affects_all <- function(df_gene, printLog, ZP_env){
  func_start(ZP_env)
  concern_info <- df_gene %>% 
    filter(vn_status==2) %>%
    .[1,] %>%
    mutate(n_mut=nrow(df_gene),
           info=paste("all copies affected by variant:", mut_id),
           conf=1,
           eval_by=ifelse(str_detect(class, "homdel"), "homdel","aff_cp")) %>%
    select(gene, n_mut, score=vn_status, conf, info, eval_by, wt_cp) %>%
    mutate(warning=NA, wt_cp_range=NA, phasing=NA)
  append_loglist(concern_info$info, ZP_env=ZP_env)
  func_end(ZP_env)
  return(concern_info)
}
#' @keywords internal
#' if no variant affecta all copies
#' @importFrom dplyr mutate select
eval_one_mut <- function(df_gene, printLog, ZP_env){
  func_start(ZP_env)
  concern_info <- df_gene %>%
    mutate(n_mut=1,
           conf=1,
           eval_by="aff_cp") %>%
    select(gene, n_mut, score=vn_status, conf, 
           info=pre_info, 
           wt_cp,
           eval_by) %>%
    mutate(warning=NA, wt_cp_range=NA, phasing=NA)
  append_loglist(concern_info$info, ZP_env=ZP_env)
  func_end(ZP_env)
  return(concern_info)
}
eval_lost_in_tumor <- function(df_gene, printLog, ZP_env){
  func_start(ZP_env)
  concern_info <- df_gene %>%
    mutate(n_mut=1,
           conf=1,
           eval_by="allele-freq",
           score=1) %>%
    select(gene, n_mut, score, conf, 
           info=pre_info, 
           wt_cp,
           eval_by) %>%
    mutate(warning="variant lost in tumor", wt_cp_range=NA, phasing=NA)
  append_loglist(concern_info$info, ZP_env=ZP_env)
  func_end(ZP_env)
  return(concern_info)
}
#' @importFrom dplyr mutate select filter left_join bind_rows
remove_duplicated_variants <- function(df_gene_raw, ZP_env){
  vm(as.character(sys.call()[1]), 1, ZP_env)
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
  func_end(ZP_env)
  return(df_gene)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom dplyr filter mutate as_tibble case_when
predict_zygosity_genewise <- function(GENE, 
                                      evaluation_per_variant, 
                                      bamDna, 
                                      bamRna, 
                                      showReadDetail,
                                      printLog, 
                                      purity, 
                                      sex, 
                                      haploBlocks,
                                      vcf, 
                                      distCutOff, 
                                      logDir, 
                                      somCna, 
                                      snpQualityCutOff, 
                                      phasingMode,
                                      AllelicImbalancePhasing, 
                                      ZP_env){
  #func_start(ZP_env)
  vm("predict_zygosity_genewise", 1, ZP_env) ## do not replcae!!!
  start_gene_eval <- Sys.time()
  increment_loglist(ZP_env)
  append_loglist(GENE, ZP_env=ZP_env)
  gene <- pre_info <-  chr <-  pos <- wt_cp <-  mut_id <- status <- . <- 
    score <- comb <- dist <- tcn <- info <- n <- all_comb <- 
    read_level_phasing_info <- copy_number_phasing_info <- mat_phased_gene <- 
  mat_info_gene <- NULL 
  #start_gene_eval <- Sys.time()
  pre_df_gene <- evaluation_per_variant %>% filter(gene==GENE)
  ## all germline variants which are lost in the tumor are excluded, 
  ## vn_status = -1
  ## check if variants at the same position are present
  df_gene_raw <- remove_duplicated_variants(pre_df_gene, ZP_env=ZP_env)
  df_gene <- df_gene_raw %>%
    filter(vn_status>=0)
 #print(10)
 #print(df_gene_raw)
 #print(df_gene)
  ## check if variants remain in that gene
  if(nrow(df_gene)==0){
    eval_for_gene <- eval_lost_in_tumor(df_gene_raw, printLog, ZP_env=ZP_env)
   #print(eval_for_gene)
      ## all variants lost in tumor
      append_loglist("all variants lost in tumor", ZP_env=ZP_env)
  } else if(any(df_gene$vn_status==2)){
      append_loglist("one of", nrow(df_gene), "variants affects all copies", ZP_env=ZP_env)
      eval_for_gene <- eval_one_mut_affects_all(df_gene, printLog, ZP_env)
    ## (2): if only one variant is present  
  } else if(nrow(df_gene)==1){
      append_loglist("one variant that does not affect all copies", ZP_env=ZP_env)
      eval_for_gene <- eval_one_mut(df_gene, printLog, ZP_env)
    ## (3): more than one variant present in gene   
  } else {  
      append_loglist(nrow(df_gene),
            "heterozygous variants detected: Initializing haplotype phasing", ZP_env=ZP_env)
      if(!is.null(logDir)){
        geneDir <- file.path(logDir, GENE)
        dir.create(geneDir)
      } else {
        geneDir <- NULL
      }
      ## returns two tibbles: first one direct phasing combinations
      ## second one indirect phasing combinations
      full_phasing_result <- phase(df_gene, somCna, bamDna, purity, sex, bamRna, 
                                   haploBlocks, vcf, distCutOff, printLog,
                                   geneDir, showReadDetail, 
                                   snpQualityCutOff, phasingMode, 
                                   AllelicImbalancePhasing, ZP_env)
     #print(full_phasing_result)
      all_comb <- full_phasing_result[[1]] %>%
        mutate(gene=GENE)
      eval_for_gene <- eval_phasing_new(all_comb, df_gene,  
                                    printLog, ZP_env)
      #print(full_phasing_result)
      read_level_phasing_info <- full_phasing_result[[2]] %>%
        mutate(gene=GENE)
      copy_number_phasing_info <- full_phasing_result[[5]] %>%
        mutate(gene=GENE)
      if(nrow(copy_number_phasing_info)==0){
        copy_number_phasing_info <- NULL
      }
      store_log(geneDir, read_level_phasing_info, "all_phasing_combinations.tsv")
      store_log(geneDir, copy_number_phasing_info, "all_phasing_combinations.tsv")
      mat_phased_gene <- list()
      mat_info_gene <- list()
      mat_phased_gene[[GENE]] <- full_phasing_result[[3]]
      mat_info_gene[[GENE]] <- full_phasing_result[[4]]
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
                               all_comb,
                               read_level_phasing_info,
                               mat_phased_gene,
                               mat_info_gene,
                               copy_number_phasing_info)
                          )
  log_message <- unlist(ZP_env$global_ZygosityPredictor_variable_debug_loglist) %>% paste(collapse = "\n")
  func_end(ZP_env)
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
                                             templateGenes, ZP_env){
  func_start(ZP_env)
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
  func_end(ZP_env)
  return(combined_uncovered)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr as_tibble group_by mutate ungroup filter
#' @importFrom purrr compact
combine_main_variant_tables <- function(df_germ, df_som, df_homdels,
                                        templateGenes, purity, ZP_env){
  func_start(ZP_env)
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
      "Please provide an annotation that contains every gene from the",
      "variant inputs"
    )
  }
  func_end(ZP_env)
  return(df_all_mutations)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr bind_rows group_by mutate relocate select ungroup relocate
bind_incdel_to_pre_eval <- function(df_incompletedels, df_all_mutations, ZP_env){
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
bind_incdel_to_final_eval <- function(df_incompletedels, final_output, ZP_env){
  func_start(ZP_env)
  gene <- 0
  if(is.null(final_output)&is.null(df_incompletedels)){
    full_output <- NULL
  } else {
    if(is.null(df_incompletedels)){
      full_output_pre <- final_output
    } else {
      full_output_pre <- bind_rows(
        final_output,
        df_incompletedels %>%
          filter(!gene %in% final_output$gene) %>%
          select(gene) %>%
          unique() %>%
          mutate(status="wt_copies_left",
                 info="somatic-incompletedel")
      )       
    }
    full_output <- full_output_pre %>%
      mutate(status=factor(status, 
                           levels=c("all_copies_affected", "wt_copies_left", 
                                    "undefined")))
  }
  func_end(ZP_env)
  return(full_output)
}

set_global_variables <- function(debug, verbose, printLog, ZP_env){
  ZP_env$global_ZygosityPredictor_variable_debug_call_depth <- 0
  ZP_env$global_ZygosityPredictor_variable_debug_log_depth <- 0
  ZP_env$global_ZygosityPredictor_variable_debug_timelist  <- list()
  ZP_env$global_ZygosityPredictor_variable_debug_verbose <- verbose
  ZP_env$global_ZygosityPredictor_variable_debug_printLog <- printLog
}
#' @keywords internal
catt <- function(printLog=FALSE, level, text){
  if(printLog==TRUE){
    message(rep("  ", level), text)
  }
}
increment_loglist <- function(ZP_env){
  ZP_env$global_ZygosityPredictor_variable_debug_loglist <- list()
}
#' @importFrom stringr str_replace_all %>%
append_loglist <- function(..., ZP_env){
  filler <- paste(rep("  ", ZP_env$global_ZygosityPredictor_variable_debug_call_depth), collapse = "")
  appendix <- paste(..., collapse = " ") %>%
    paste0(filler,.) %>%
    str_replace_all("\n", paste0("\n",filler))
  ZP_env$global_ZygosityPredictor_variable_debug_loglist <- append(ZP_env$global_ZygosityPredictor_variable_debug_loglist, appendix)
  if(ZP_env$global_ZygosityPredictor_variable_debug_printLog==TRUE){
    message(appendix)
  }
}
add_timestamp <- function(ZP_env){
  ZP_env$global_ZygosityPredictor_variable_debug_timelist[[ZP_env$global_ZygosityPredictor_variable_debug_call_depth]] <- Sys.time()
}
remove_timestamp <- function(ZP_env){
  ZP_env$global_ZygosityPredictor_variable_debug_timelist <- 
    ZP_env$global_ZygosityPredictor_variable_debug_timelist[c(1:ZP_env$global_ZygosityPredictor_variable_debug_call_depth)]
}
#' @importFrom magrittr %>%
vm <- function(mes, depth=0, ZP_env){
  if(depth>0){
    to_add <- "|-"
  }  
  if(depth<0){
    to_add <- "+-"
  }
  if(ZP_env$global_ZygosityPredictor_variable_debug_verbose==TRUE){
    to_print <- paste0(paste(rep("| ", ZP_env$global_ZygosityPredictor_variable_debug_call_depth),collapse=""), to_add,mes)
  } 
  if(depth>0){
    increment_call_depth(ZP_env)
    add_timestamp(ZP_env)
    timestamp <- ""
  }  
  if(depth<0){
    timestamp <- pt(ZP_env$global_ZygosityPredictor_variable_debug_timelist[[ZP_env$global_ZygosityPredictor_variable_debug_call_depth]], 
                    Sys.time()) %>%
      paste("  ~",.)
    decrement_call_depth(ZP_env)
    remove_timestamp(ZP_env)
  }
  if(ZP_env$global_ZygosityPredictor_variable_debug_verbose==TRUE){
    message(to_print, timestamp)
  }
}