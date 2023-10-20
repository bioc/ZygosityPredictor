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
#' @importFrom stringr %>% str_split
#' @importFrom dplyr bind_rows mutate_at
#' @importFrom purrr compact
phase_along_path <- function(path, sub_checked_read_presence, bam_raw,
                             purity, printLog){
  . <- mut_id1 <- mut_id2 <- comb_id <- NULL
  rel_sub_combs <- lapply(seq_len(length(path)), function(N){
    if(N==length(path)){
      return(NULL)
    } else {
      sub_checked_read_presence %>%
        filter(
          (path[[N]]==mut_id1|path[[N]]==mut_id2)&
            (path[[N+1]]==mut_id1|path[[N+1]]==mut_id2)
        ) %>%
        return()
    }
  }) %>%
    compact() %>%
    bind_rows() %>%
    mutate_at(.vars = c("chr1", "chr2", "pos1", "pos2", "distance"),
              .funs=as.numeric)
  bam <- bam_raw %>% compact() %>% Reduce(function(x,y)c(x,y),.)
  sub_fin_comb <-  apply(rel_sub_combs, 1, function(sub_comb){
    catt(printLog, 6, c("sub combination:", sub_comb[["comb_id"]], 
                        "distance:", sub_comb[["distance"]]))
    bam_fil <- bam %>% .[which(.$comb_id == sub_comb[["comb_id"]])]
    sub_classified_reads <- lapply(unique(bam_fil$qname), 
                                   core_tool,
                                   bam_fil,
                                   as.numeric(sub_comb[["pos1"]]),
                                   as.numeric(sub_comb[["pos2"]]),
                                   sub_comb[["alt1"]],
                                   sub_comb[["alt2"]],
                                   sub_comb[["ref1"]],
                                   sub_comb[["ref2"]],
                                   sub_comb[["class1"]],
                                   sub_comb[["class2"]]) %>%
      compact() %>%
      bind_rows()
    sub_status_combination <- 
      classify_combination(sub_classified_reads,
                           purity,
                           eval_full=FALSE,
                           printLog)
    return(sub_status_combination %>%
             mutate(comb_id=sub_comb[["comb_id"]]))
  }) %>% bind_rows()%>%
    mutate(
      f=str_split(comb_id, "-") %>% map_chr(.,1),
      l=str_split(comb_id, "-") %>% map_chr(.,2),
    )
  return(sub_fin_comb)
}

perform_extended_snp_phasing <- function(df_gene, vcf, refGen, verbose){
  
  region_to_load <- paste0(unique(df_gene$chr), ":", min(df_gene$pos),
                           "-", max(df_gene$pos)) %>%
    GRanges()
  print(paste0(unique(df_gene$chr), ":", min(df_gene$pos),
               "-", max(df_gene$pos)))
  loaded_vcf <- loadVcf(vcf, unique(df_gene$chr), region_to_load, refGen, 
                        verbose)
  print(loaded_vcf)
}


perform_extended_snp_phasing_old <- function(main_comb, vcf,
                                         df_gene, bamDna, bamRna, 
                                         purity, full_read_info, RESULT,
                                         ext_snp_phasing,
                                         printLog, verbose){
  catt(printLog, 4, "Initializing extended SNP phasing")
  vars_in_between_raw <- check_for_snps_between_main_muts(main_comb, vcf,
                                                          df_gene, 
                                                          snp_dist=0) %>%
    filter(site=="between") %>%
    select(-site)
  
  if(nrow(vars_in_between_raw)>2){
    catt(printLog, 5, c(nrow(vars_in_between_raw)-2, 
                        "SNPs between variants detected"))
    vars_in_between <- vars_in_between_raw %>% 
      mutate(mut_id=
               c(mut_id[which(!is.na(mut_id))],
                 paste0("s", 
                        seq_len(nrow(.)-length(which(!is.na(mut_id)))))))
    sub_combinations <- make_phasing_combinations(vars_in_between) %>%
      .[which(!(str_detect(.$mut_id1, "m")&str_detect(.$mut_id2, "m"))),]
    sub_checked_read_presence_list <- apply(sub_combinations, 1, 
                                            check_read_presence, 
                                            bamDna, bamRna)
    sub_checked_read_presence <- lapply(sub_checked_read_presence_list,
                                        nth, n=1) %>%
      bind_rows() %>%
      filter(nreads>0)
    bam_raw <- lapply(sub_checked_read_presence_list,
                      nth, n=2)
    ## now check with the remaining ones if there is a connection
    ## first check if the main variants are still there
    still_present_muts <- c(sub_checked_read_presence$mut_id1,
                            sub_checked_read_presence$mut_id2) %>%
      unique()
    if(main_comb[["mut_id1"]] %in% still_present_muts&
       main_comb[["mut_id2"]] %in% still_present_muts){
      
      ## make list to use with A* algorithm
      
      #print(still_present_muts)
      #print(sub_checked_read_presence)
      #print(main_comb)
      
      # path <- find_path(still_present_muts, sub_checked_read_presence,
      #                   main_comb)
      path <- find_path(sub_checked_read_presence %>%
                          select(mut_id1, mut_id2), main_comb[["mut_id1"]],
                        main_comb[["mut_id2"]])
      if(!is.null(path)){
        catt(printLog, 5, c("path found:", 
                            paste(unlist(path), collapse=" - ")))
        ## find relevant combinations to phase
        if(printLog==TRUE){
          apply(
            vars_in_between[
              which(vars_in_between$mut_id %in% unlist(path)),],
            1, function(r){
              catt(printLog, 5, c(r[["mut_id"]],
                                  r[["chr"]],
                                  r[["pos"]]))
            })
        }
        sub_fin_comb <- phase_along_path(path, sub_checked_read_presence,
                                         bam_raw, purity, printLog)
        non_null <- sub_fin_comb %>%
          filter(status!="null")
        if(nrow(non_null)!=0&
           (main_comb[["mut_id1"]] %in% c(non_null$f, non_null$l)&
            main_comb[["mut_id2"]] %in% c(non_null$f, non_null$l))){
          catt(printLog, 5, 
               "all combinations of path have assigned status")
          final <- recombine_phasing_results(non_null, main_comb)
          ext_snp_phasing <- sub_fin_comb %>%
            mutate(main_comb=main_comb[['comb_id']])
          if(nrow(final)!=0){
            return_null_result <- FALSE
            RESULT <- finalize_snp_phasing(final, main_comb)
            full_read_info <- NULL    
            ext_snp_info <- paste(
              "SNP phasing successfull:", RESULT$status, "; wt cp:",
              RESULT$wt_cp)
          } else {
            ext_snp_info <- 
              "recombination impossible: more than one diff status in path"
          }
        } else {
          ext_snp_info <- "recombination imposiible: null status in path"
        }
      } else {
        ext_snp_info <- "no connecting path between main muts"
      }
    } else {
      ext_snp_info <- "no overlap with main muts"
    }
  } else {
    ext_snp_info <- "no SNPs between"
  }
  return(list(RESULT=RESULT, full_read_info=full_read_info, 
              ext_snp_info=ext_snp_info, ext_snp_phasing=ext_snp_phasing))
}


#' @keywords internal
#' @importFrom stringr %>% str_count str_replace_all str_detect
#' @importFrom dplyr desc
recombine_phasing_results <- function(non_null, main_comb){
  status <- comb_id <- f <- l <- fin_comb <- NULL
  sames <- 
    non_null %>%
    arrange(desc(status))  %>%
    select(status, comb_id, f, l)
  new <- sames
  for(n in seq_len(str_count(sames$status, "same") %>% sum())){
    curr_comb_id <- new[[n, "comb_id"]]
    to_rep <- paste(paste0(new[[n, "f"]],"$"),
                    paste0(new[[n, "l"]],"$"),
                    sep="|")
    replacement <- paste0(new[[n, "f"]],new[[n, "l"]])
    new <- new %>% 
      mutate(f=str_replace_all(f, to_rep, replacement),
             l=str_replace_all(l, to_rep, replacement)) 
  }
  final <- new %>%
    mutate(fin_comb=paste(f, l, sep="-")) %>%
    filter(str_detect(fin_comb, main_comb[["mut_id1"]])&
             str_detect(fin_comb, main_comb[["mut_id2"]]))
  return(final)
}
#' @keywords internal
#' @importFrom dplyr tibble
finalize_snp_phasing <- function(final, main_comb){
  status <- unique(final$status)
  tcn1 <- main_comb[["tcn1"]]
  tcn2 <- main_comb[["tcn2"]]
  aff_copies1 <- main_comb[["aff_cp1"]]
  aff_copies2 <- main_comb[["aff_cp2"]]
  psc <- pre_scoring(tcn1, tcn2, status, aff_copies1, aff_copies2)
  RESULT <- tibble(
    comb=main_comb[['comb_id']],
    class_comb=paste(main_comb[['class1']], 
                     main_comb[['class2']], sep="-"),
    dist=as.numeric(main_comb[["distance"]]),
    score=psc$pre_score,
    wt_cp=round(psc$left_wt_copies,2),
    tcn=round(psc$mtcn,2),
    info="status defined after extended HP phasing",
    status=status,
    DNA_rds=0,
    RNA_rds=0,
    both=0,
    none=0,
    mut1=0,
    mut2=0,
    no_overlap=0,
    none_raw=0,
    dev_var=0
  )
  return(RESULT)
}

