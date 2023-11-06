phase_combination <- function(mat_gene_relcomb, comb, bamDna, bamRna, verbose, geneDir, 
                              phasingDir, phasing_type, showReadDetail){
  
  classified_main_comb <- tibble(comb=comb,
                                 status="null",
                                 nstatus=0,
                                 conf=0,
                                 phasing=comb)
  main_classified_reads <- classify_reads(mat_gene_relcomb, bamDna, bamRna, verbose)
  append_loglist(nrow(main_classified_reads), 
                 "reads / read-pairs covering both positions")
  if(nrow(main_classified_reads)!=0){
    #return_null_result <- FALSE
    if(showReadDetail==TRUE){
      store_log(geneDir, main_classified_reads %>%
                  mutate(comb=comb), 
                paste0("classified_reads_", comb,".tsv"))
    } 
    classified_main_comb <- classify_combination(main_classified_reads,
                                                 purity,
                                                 printLog,
                                                 verbose
    ) %>%  
      mutate(
        dist=mat_dist[comb],
        class_comb=paste(mat_gene_relcomb[,"class"][[2]], mat_gene_relcomb[,"class"][[1]], 
                         sep="-"),
        comb=comb,
        phasing=comb
      ) 
    
    store_log(phasingDir, classified_main_comb %>% mutate(gene=basename(geneDir)), 
              paste0("classified_combination_", 
                     comb, ".tsv"))
  }   
  append_matrices(classified_main_comb)
  return(classified_main_comb)
}


perform_direct_phasing <- function(mat_gene, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail, geneDir){
  func_start()
  phasingDir <- NULL
  phasing_type <- "direct"
  if(!is.null(geneDir)){
    phasingDir <- file.path(geneDir, phasing_type)
    dir.create(phasingDir)
  }
  
  phasing_info <- tibble()
  unphased <- which(is.na(mat_phased))
  while(length(unphased)>0){
    #print(unphased)
    ## as long as there are unphased combinations, try to phase
    relcombxy <- get_xy_index(unphased[1], nrow(mat_phased))
    mat_gene_relcomb <- mat_gene[relcombxy,]
    classified_main_comb <- phase_combination(mat_gene_relcomb, unphased[1], bamDna, bamRna, 
                                              verbose, geneDir, phasingDir, phasing_type, showReadDetail)
    phasing_info <- bind_rows(phasing_info,
                              classified_main_comb)
    unphased <- which(is.na(mat_phased))
  }
  
  
  
  
  
  
  # if(nrow(phasing_info)>0){
  #   ## if something could be phased, check if unphased combinations can be
  #   ## solved by master_combinations
  #   main_master_combs <- subdivide_muts_in_three_combs(phasing_info$comb)
  #   
  #   if(!is.null(main_master_combs)){
  #     phasing_status_filled <- phasing_info %>%
  #       select(comb, nstatus, conf, phasing) %>%
  #       define_next_priority(main_master_combs, .) %>%
  #       iterate_master_combs() %>%
  #       ungroup() %>%
  #       select(-master_comb) %>% unique()
  #   } else {
  #     phasing_status_filled <- phasing_info 
  #   }
  #   phasing_status <- phasing_status_filled %>%
  #     mutate(nstatus=ifelse(is.na(nstatus), 0, nstatus),
  #            conf=ifelse(is.na(conf), 0, conf),
  #            status=get_string_status(nstatus)) %>%
  #     select(comb, status, nstatus, conf, phasing) 
  # } else {
  #   phasing_status <- NULL
  # }
  # filled_phasing_status <- fill_phasing_status(all_combinations, phasing_status, 
  #                                              "direct", verbose)
  func_end()
  return(list(status=NULL,
              info=phasing_info,
              exit=tibble(phasing="direct", info="test")))
}
