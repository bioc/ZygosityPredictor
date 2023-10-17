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
      vm("defining result for combination", verbose)
      RESULT <- classify_combination(main_classified_reads,
                                     purity,
                                     eval_full = TRUE, 
                                     printLog,
                                     main_comb[['tcn1']],
                                     main_comb[['tcn2']],
                                     main_comb[['aff_cp1']],
                                     main_comb[['aff_cp2']],
                                     verbose
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
      vm("result is null", verbose)
      RESULT <- return_null_result(main_comb, bamRna, "test")
      full_read_info <- NULL
    }
    vm("returning result", verbose)
    return(list(assigned=RESULT, detailed=full_read_info, 
                #snp_phasing=ext_snp_phasing, 
                level_2_phasing=level_2_phasing))
  })
  return(final_combinations)
}