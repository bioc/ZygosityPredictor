perform_direct_phasing <- function(all_combinations, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail, geneDir){
  func_start(sys.call(), verbose)
  phasing_info <-  apply(all_combinations, 1, function(main_comb){
    catt(printLog, 2 ,c("combination:", main_comb[["comb_id"]]))
    #return_null_result <- TRUE
    RESULT <- full_read_info <- NULL
    main_classified_reads <- classify_reads(main_comb, bamDna, bamRna, verbose)
    catt(printLog, 3, 
         c(nrow(main_classified_reads), 
           "reads / read-pairs overlapping both positions"))
    if(nrow(main_classified_reads)!=0){
      #return_null_result <- FALSE
      if(showReadDetail==TRUE){
        store_log(geneDir, main_classified_reads %>%
          mutate(comb=main_comb[['comb_id']]), 
          paste0("classified_reads_",main_comb[['comb_id']],".tsv"))
      } 
      RESULT <- classify_combination(main_classified_reads,
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
      if(RESULT$status!="null"){
        #return_null_result <- FALSE
        catt(printLog, 3, c("status:", RESULT$status))
      } else {
        catt(printLog, 3, "status can not be defined")
      }
    }   
    # if(return_null_result==TRUE){
    #   RESULT <- return_null_result(main_comb, bamRna, "test")
    # }
    return(RESULT)
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
              exit="test"))
}