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
  func_start()
  graph <- graph_from_data_frame(connections, directed = FALSE)
  shortest_path <- 
    all_shortest_paths(graph, from = mut_id1, to = mut_id2, mode = "all")$res
  if (length(shortest_path) > 0) {
    path_list <- lapply(shortest_path, names)
  } else {
    path_list <- NULL
  }
  return(path_list)
  func_end()
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

get_muts_to_be_phased <- function(df_gene){
  warning("IMPLEMENT FILTER OF RELEVANT PHASING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
  return(df_gene$mut_id)
}


rearrange_priorities <- function(sub_somb, phasing_status){
  
  prioritized <- 
    sub_comb %>%
    left_join(phasing_status, by=c("comb_id"="comb")) %>%
    mutate(isna=is.na(nstatus)) %>%
    arrange(desc(isna))
  
  return(prioritized)
  
}

load_all_snps_in_relevant_gene_region <- function(df_gene, phasedVcf, refGen, verbose, distCutOff){
  func_start()
  region_to_load <- paste0(unique(df_gene$chr), ":", min(df_gene$pos)-distCutOff,
                           "-", max(df_gene$pos)+distCutOff) %>%
    GRanges()  %>%
    split_genomic_range(.,df_gene$pos)
  
  loaded_vcf_hz <- loadVcf(phasedVcf, unique(df_gene$chr), region_to_load, refGen, 
                           verbose, "HZ") 
  #define_region_to_load(df_gene)
  snps <- loaded_vcf_hz %>% as_tibble() %>%
    mutate(snp_id=paste0("s", c(1:nrow(.))))
  
  append_loglist(nrow(snps),"SNPs detected between main mutations")
  sub_comb <- combine_main_muts_and_snps(df_gene, snps) %>% 
    #filter(dist<distCutOff) %>%
    arrange(mut_id1)
  func_end()
  return(sub_comb)
}

recombine_muts <- function(current_path, M){
  outer(sort(current_path), sort(current_path), `paste`) %>% 
    .[which(upper.tri(.))] %>% str_replace_all(" ", "-") %>%
    return()
  # lapply(c(2:length(current_path)), function(i){
  #   paste(sort(c(current_path[i-1], current_path[i])), collapse="-")
  # }) %>% unlist() %>% c(.,M) %>%
  #   return()
}

append_main_pashing <- function(phasing_result_df, classified_main_comb){
  func_start()
  new_result <- classified_main_comb %>% select(comb, nsnew=nstatus, 
                                                cnew=conf)
  new_phasing_result_df <- 
    left_join(
    phasing_result_df %>% select(comb, nsold=nstatus, cold=conf),
    new_result,
    by="comb"
  ) %>% rowwise() %>%
    mutate(nstatus=case_when(
      is.na(nsold)&is.na(nsnew) ~ NA,
      is.na(nsold) ~ nsnew,
      is.na(nsnew) ~ nsold,
      TRUE ~ max(nsnew, nsold)
    ),
    conf=case_when(
      is.na(nsold)&is.na(nsnew) ~ NA,
      is.na(nsold) ~ cnew,
      is.na(nsnew) ~ cold,
      nsnew==nsold ~ max(cnew, cold),
      nsnew>nsold ~ cnew,
      TRUE ~ cold
    )) %>%
      select(comb, nstatus, conf)
  func_end()
  return(new_phasing_result_df)
}

iterate_master_combs <- function(new_res){
  func_start()
  inp <- new_res %>% group_by(master_comb)
  something_changed <- TRUE
  while(something_changed==TRUE){
   #print(1)
    out <- inp %>%
      mutate(nstatus=recombine_three_combinations(nstatus),
             conf=case_when(
               (is.na(conf)|conf==0)&sum(conf==0|is.na(conf))==1 ~ mean(conf[which(conf!=0)]),
               TRUE ~ conf))
    inp_na_rm <- inp$nstatus
    out_na_rm <- out$nstatus
    inp_na_rm[which(is.na(inp_na_rm))] <- 0
    out_na_rm[which(is.na(out_na_rm))] <- 0
    
    if(sum(inp_na_rm==out_na_rm)==length(inp_na_rm)){
      something_changed <- FALSE
     #print(out)
    } else {
      inp <- inp %>%
        select(-nstatus, -conf) %>%
        left_join(
          ungroup(out) %>% select(comb, nstatus, conf) %>% unique() %>% filter(!(is.na(nstatus)|nstatus==0)),
          by="comb"
        )
     #print(inp)
    }
  }
  func_end()
  return(out)
  
}
next_path <- function(sub_comb, phasing_result_df, muts){
  null_combs <- phasing_result_df[which(phasing_result_df$nstatus==0),] %>%
    pull(comb)
  find_path(sub_comb %>%
              filter(!comb_id %in% null_combs) %>%
              select(mut_id1, mut_id2), muts[1], muts[2]) %>%
    return()
}
finalize_indirect_phasing <- function(all_combinations_to_be_phased, states_main_comb){
  func_start()
  all_combinations_to_be_phased %>%
    filter(!any(nstatus==0))
  
  final <- states_main_comb %>%
    left_join(all_combinations_to_be_phased %>% group_by(comb) %>% summarize(conf=max(conf))) %>%
    mutate(phasing=phasing_type,
           status=case_when(nstatus==2 ~ "diff",
                            nstatus==1 ~ "same",
                            TRUE ~ "null"
           ))
  func_end()
  return(final)
}
perform_indirect_phasing <- function(df_gene, direct_phasing, phasedVcf, refGen, 
                                     verbose, phasingMode, all_combinations, showReadDetail, geneDir, distCutOff){
  
  func_start()
  phasing_type <- "indirect"
  do_indirect_phasing <- decide_following_phasing(list(direct_phasing), 
                                                  phasedVcf, verbose,
                                       phasingMode, phasing_type)
  
  
  if(do_indirect_phasing){
    append_loglist("Initialize indirect phasing")
    if(!is.null(geneDir)){
      phasingDir <- file.path(geneDir, phasing_type)
      dir.create(phasingDir)
    }
    main_comb_to_be_phased <- direct_phasing$status %>%
      filter(nstatus==0|is.na(nstatus)) %>%
      pull(comb)
    append_loglist("Initialize indirect phasing of combinations:\n", 
                   paste(main_comb_to_be_phased, collapse = ", "))
    sub_comb_raw <- load_all_snps_in_relevant_gene_region(df_gene, phasedVcf, 
                                                      refGen, verbose, 
                                                      distCutOff=0)
    sub_comb <- sub_comb_raw %>% filter(dist<distCutOff)

    present_main_muts <- c(sub_comb$mut_id1, sub_comb$mut_id2) %>%
      .[which(str_detect(.,"m"))] %>% unique()
    ## safe input phasing results
    input_phasing_results <- direct_phasing$status %>% select(comb, nstatus, conf)
    main_master_combs <- subdivide_muts_in_three_combs(input_phasing_results$comb)
    ## define global overview about combinations that have been phased already
    phasing_result_df <- input_phasing_results %>%
      bind_rows(
        sub_comb_raw %>%
          filter(!comb_id %in% input_phasing_results$comb) %>%
          mutate(nstatus=ifelse(dist<=distCutOff, NA, 0),
                 conf=ifelse(dist<=distCutOff, NA, 0)) %>%
          select(comb=comb_id, nstatus, conf))
    ## now check by main combination
    MCC <- 1
    all_main_combs_solved <- FALSE
    while(all_main_combs_solved==FALSE&MCC<=length(main_comb_to_be_phased)){
      ## get muts from main combination
      M <- main_comb_to_be_phased[MCC]
      append_loglist("Main phasing of", M)
      muts <- unlist(str_split(M, "-"))
      ## check if they are both present in any of the combinations
      if(length(intersect(muts, present_main_muts))==2){
        all_paths <- next_path(sub_comb, phasing_result_df, muts)
        if(!is.null(all_paths)){
          unchecked_paths <- length(all_paths)
          main_mut_reached <- FALSE
          combined_phasing_results <- tibble()
          while(unchecked_paths>0&main_mut_reached==FALSE){
            append_loglist("possible path detected:", 
                           paste(all_paths[[unchecked_paths]],collapse="->"))
            ## calculate master combinations for this path
            master_combs <- subdivide_muts_in_three_combs(recombine_muts(all_paths[[unchecked_paths]]))
            ## calculate phasing priorities for which should be the next 
            ## combination to be phased
            all_combinations_to_be_phased <- define_next_priority(master_combs, phasing_result_df)
            max_combs_to_phase <- all_combinations_to_be_phased %>%
              filter(is.na(nstatus)) %>% nrow()
            ## predifine variables for phasing loop
            to_phase_next <- pick_next(all_combinations_to_be_phased) 
            i <- 1
            while(main_mut_reached==FALSE&i<=max_combs_to_phase&sum(all_combinations_to_be_phased$prio)>0){
              append_loglist("phasing of", to_phase_next)
              ## do phasing of prio 1 combination
              classified_main_comb <- phase_combination(
                sub_comb[which(sub_comb$comb_id==to_phase_next),], bamDna, bamRna, 
                verbose, geneDir, phasingDir, phasing_type, showReadDetail)
              ## store phasing info for combination
              combined_phasing_results <- bind_rows(combined_phasing_results, classified_main_comb)
              ## append new result to global phasing info
              phasing_result_df <- append_main_pashing(phasing_result_df, 
                                                       classified_main_comb)
              ##  define next prioroity
              all_combinations_to_be_phased <- define_next_priority(master_combs, phasing_result_df) %>%
                ## and try to solve master combinations
                iterate_master_combs()
              ## append global phasing results (if something could be solved by master_comb)
              phasing_result_df <- append_main_pashing(phasing_result_df, 
                                                       unique(ungroup(all_combinations_to_be_phased[c(1,4,5)])))
              append_loglist(print_tibble(phasing_result_df))
              ## check if main combination could be solved
              states_main_comb <- all_combinations_to_be_phased %>%
                filter(comb %in% main_comb_to_be_phased) %>%
                group_by(comb) %>%
                summarize(nstatus=max(nstatus))
              if(0 %in% states_main_comb$nstatus){
                to_phase_next <- pick_next(all_combinations_to_be_phased) 
              } else {
                main_mut_reached <- TRUE
              }
              i <- i+1
            } ## while loop over single combinations
            if(main_mut_reached==TRUE){
              append_loglist("successfully phased along path")
              ## check whioch three-comb was the relevant one
              #final <- finalize_indirect_phasing(all_combinations_to_be_phased, states_main_comb){
              exit <- NULL
            } else {
              unchecked_paths <- unchecked_paths-1
              if(unchecked_paths==0){
                ## check if another path can be found
                all_paths <- next_path(sub_comb, phasing_result_df, muts)
                if(!is.null(all_paths)){
                  unchecked_paths <- length(all_paths)
                }
              }
              ## phasing did not work       
              exit <- "phasing did not work" 
            }  
          } ## while loop over paths
          append_loglist("all paths checked")
        } else {
          exit <- "no path between main muts found" 
          append_loglist(exit)
        }
      } else {
        exit <- "at least one mut of combinations has no SNP closer than dist limit"
        append_loglist(exit)
      } 
      output_main_combs <- define_main_comb_output(main_master_combs, phasing_result_df, input_phasing_results, phasing_type)
      if(!0 %in% output_main_combs$nstatus){
        all_main_combs_solved <- TRUE
      }
      MCC <- MCC+1
    } ## for loop over main combinations
  } 
  return(output_main_combs)    
    

  #   } else {
  #     ## no path found inside distCutOff
  #   }
  # } else {
  #   final <- fill_phasing_status(all_combinations, NULL, 
  #                                phasing_type, verbose)
  # }
  # indirect_phasing_results <- list(
  #   status=final,
  #   info=combined_phasing_results,
  #   exit=NA
  # )
  # func_end()
  # return(indirect_phasing_results)
  ## aggregate confidence from used phasing
}

define_main_comb_output <- function(main_master_combs, phasing_result_df, input_phasing_results, phasing_type){
  func_start()
  if(is.null(main_master_combs)){
    output_main_combs <- phasing_result_df %>%
      filter(comb %in% input_phasing_results$comb)
  } else {
    output_main_combs <- define_next_priority(main_master_combs, phasing_result_df %>%
                                                filter(comb %in% input_phasing_results$comb)) %>%
      iterate_master_combs()        
  }
  fin <- output_main_combs %>% mutate(phasing=phasing_type, 
                                      status=get_string_status(nstatus))
  func_end()
  return(fin)
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

