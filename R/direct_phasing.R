perform_direct_phasing <- function(all_combinations, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail){
  vm(as.character(sys.call()[1]), verbose, 1)
  phasing_info <-  apply(all_combinations, 1, function(main_comb){
    catt(printLog, 2 ,c("combination:", main_comb[["comb_id"]]))
    return_null_result <- TRUE
    RESULT <- NULL
    full_read_info <- NULL
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
                                     #eval_full = TRUE, 
                                     printLog,
                                     # main_comb[['tcn1']],
                                     # main_comb[['tcn2']],
                                     # main_comb[['aff_cp1']],
                                     # main_comb[['aff_cp2']],
                                     verbose
      ) %>%  
        mutate(
          # DNA_rds=nrow(main_classified_reads %>% filter(origin=="DNA")),
          # RNA_rds=nrow(main_classified_reads %>% filter(origin=="RNA")),
          dist=as.numeric(main_comb[["dist"]]),
          class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                           sep="-"),
          comb=main_comb[['comb_id']]
        ) 
      if(RESULT$status!="null"){
        return_null_result <- FALSE
        catt(printLog, 3, c("status:", RESULT$status))
      } else {
        catt(printLog, 3, "status can not be defined")
      }
    }   
    if(return_null_result==TRUE){
      RESULT <- return_null_result(main_comb, bamRna, "test")
      full_read_info <- NULL
    }
    return(RESULT)
  }) %>%
    bind_rows() %>%
    mutate(phasing="direct")
  phasing_status <- phasing_info %>%
    select(comb, status, nstatus, conf, phasing)
  vm("  - done", verbose, -1)
  return(list(status=phasing_status,
              info=phasing_info))
}