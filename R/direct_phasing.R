perform_direct_phasing <- function(all_combinations, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail, geneDir){
  func_start(verbose)
  phasing_info <-  apply(all_combinations, 1, function(main_comb){
    append_loglist("direct phasing of combination:", main_comb[["comb_id"]],
                   "distance:", main_comb[["dist"]])
    classified_main_comb <- full_read_info <- NULL
    main_classified_reads <- classify_reads(main_comb, bamDna, bamRna, verbose)
    append_loglist(nrow(main_classified_reads), 
           "reads / read-pairs covering both positions")
    if(nrow(main_classified_reads)!=0){
      #return_null_result <- FALSE
      if(showReadDetail==TRUE){
        store_log(geneDir, main_classified_reads %>%
          mutate(comb=main_comb[['comb_id']]), 
          paste0("classified_reads_",main_comb[['comb_id']],".tsv"))
      } 
      classified_main_comb <- classify_combination(main_classified_reads,
                                     purity,
                                     printLog,
                                     verbose
      ) %>%  
        mutate(
          dist=as.numeric(main_comb[["dist"]]),
          class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                           sep="-"),
          comb=main_comb[['comb_id']],
          phasing="direct"
        ) 
    }   
    return(classified_main_comb)
  }) %>%
    compact() %>%
    bind_rows() 
  if(nrow(phasing_info)>0){
    phasing_status <- phasing_info %>%
      select(comb, status, nstatus, conf, phasing)    
  } else {
    phasing_status <- NULL
  }
  filled_phasing_status <- fill_phasing_status(all_combinations, phasing_status, 
                                               "direct", verbose)
  func_end(verbose)
  return(list(status=filled_phasing_status,
              info=phasing_info,
              exit=tibble(phasing="direct", info="test")))
}