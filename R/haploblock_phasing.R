phase_per_variant <- function(df_dist, df_gene, verbose, bamDna, bamRna, 
                              purity,
                              printLog){
  phased_per_mut_full <- lapply(unique(df_dist$mut_id), function(MUT_ID){
    df_dist_mut <- df_dist %>%
                     filter(mut_id==MUT_ID)
    high_conf_result <- FALSE
    i <- 1
    phased_muts <- tibble()
    ## only for documentation
    sub_status_list <- tibble()
    while(high_conf_result==FALSE&i<=nrow(df_dist_mut)){
      comb <- df_dist_mut[i,]
      sub_comb <- bind_rows(
                    df_gene %>% 
                      filter(mut_id==comb$mut_id) %>%
                      select(chr, pos, mut_id, ref, alt, class),
                    comb %>%
                      mutate(ALT=unlist(
                        lapply(ALT, function(x){as.character(x[[1]][1])})),
                          class="snv") %>%
                          select(chr=seqnames, pos=start, mut_id=snp_id, 
                                 ref=REF, alt=ALT, class)
                        ) %>%
        make_phasing_combinations()
 ## the function classify_reads is usually used inside of apply()
## thats why the input needs to be adjusted in the format of apply iterations:
      comb_vec <- as.character(sub_comb) %>% set_names(nm=names(sub_comb))
      res <- classify_reads(comb_vec, bamDna, bamRna, verbose)
      if(nrow(res)!=0){
        sub_status_combination <- classify_combination(res,
                                                       purity,
                                                       printLog,
                                                       verbose) %>%
          mutate(dist=comb$dist,
                 class_comb=
                   paste(df_gene[which(df_gene$mut_id==comb$mut_id),]$class,
                         "snp",
                         sep="-"))
        sub_status_list <- bind_rows(sub_status_list,
                                     sub_status_combination %>% 
                                       mutate(comb=comb$comb_id))
        if(!sub_status_combination$status=="null"){
          phased <- tibble(mut_id=comb$mut_id,
                           gt=get_genotype(comb$gt, 
                             sub_status_combination$status))
          phased_muts <- bind_rows(phased_muts, phased)
          high_conf_result <- TRUE
        } 
      }
      i <- i +1 
    } 
    return(list(phased_muts, sub_status_list)) 
  }) 
  return(phased_per_mut_full)
}
annotate_haploblocks_to_variants <- function(df_gene, haploBlocks, verbose){
  ## annotate to df_gene in which haploblock the variant is located, 0 means no 
  ## haploblock
  df_gene_hap <- apply(df_gene, 1, function(mut){
    gr_roi <- GenomicRanges::GRanges(paste0(mut[["chr"]], ":", 
                                            as.numeric(mut[["pos"]]), "-",
                                            as.numeric(mut[["pos"]])))
    hap <- IRanges::subsetByOverlaps(haploBlocks, gr_roi)
    if(length(hap)==0){
      det_hap_id <- 0
    } else {
      det_hap_id <- hap$hap_id
    }
    return(as_tibble(t(mut)) %>% mutate(hap_id=det_hap_id))
  }) %>% bind_rows()
  return(df_gene_hap)
}
define_status_for_main_combination <- function(phased_per_mut, log_phasing_results){
  all_mut_ids <- phased_per_mut$mut_id
  max_conf_per_mut <- lapply(all_mut_ids, function(M){
    log_phasing_results %>%
      mutate(mut_id=str_match(comb, "m\\d+") %>% 
               unlist() %>% as.character()) %>%
      filter(mut_id==M) %>%
      filter(conf==max(conf)) %>%
      select(mut_id, conf)
  }) %>% bind_rows() %>%
    left_join(phased_per_mut, by=c("mut_id"="mut_id"))
  status_per_comb <- outer(all_mut_ids, all_mut_ids, `paste`) %>%
    .[which(upper.tri(.))] %>%
    as.data.frame() %>% 
    set_names(nm='raw') %>%
    mutate(
      mut_id1=str_split(raw, " ") %>% map_chr(.,1),
      mut_id2=str_split(raw, " ") %>% map_chr(.,2)
    ) %>%
    left_join(max_conf_per_mut %>% select(mut_id1=mut_id, gt1=gt, conf1=conf), 
              by="mut_id1") %>%
    left_join(max_conf_per_mut %>% select(mut_id2=mut_id, gt2=gt, conf2=conf), 
              by="mut_id2") %>%
    rowwise() %>%
    mutate(status=case_when(gt1==gt2 ~ "same",
                            TRUE ~ "diff"),
           nstatus=case_when(gt1==gt2 ~ 1,
                             TRUE ~ 2),
           conf=0.5*(conf1+conf2),
           comb=paste(paste0("m",
                                sort(as.numeric(str_match(
                                  c(mut_id1, mut_id2), "\\d+")))),
                         collapse = "-")
           
    ) %>%
    ungroup()
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
perform_haploblock_phasing <- function(direct_phasing, df_gene, 
                                       haploBlocks, phasedVcf,
                                    bamDna, bamRna, 
                                    purity, distCutOff,
                                    printLog, verbose,
                                    geneDir, refGen){
  per_haploblock <- log_phasing_results <- NULL
  do_haploblock_phasing <- eval_direct_results(direct_phasing, phasedVcf,
                                               haploBlocks)
  if(do_haploblock_phasing){
    catt(printLog, 4, "Initializing level 2 phasing")
    df_gene_hap <- annotate_haploblocks_to_variants(df_gene, haploBlocks, 
                                                    verbose)
    ## check if at least two variants are in the same haploblock
    if(any(df_gene_hap$hap_id > 0)&any(table(
      df_gene_hap[which(df_gene_hap$hap_id>0),]$hap_id
    )>1)){
      #vm("at least two variants detected in same haploblock", verbose, 1)
      relevant_haploblocks <- table(df_gene_hap$hap_id) %>% 
        .[which(.>1)] %>% names() %>% .[which(!.=="0")]
      ## iterate over haploblocks
      per_haploblock <- lapply(relevant_haploblocks, function(HAP_ID){
        log_phasing_results <- NULL
        region_to_load <- haploBlocks[which(haploBlocks$hap_id==HAP_ID)]
        loaded_vcf <-  loadVcf(phasedVcf, unique(df_gene$chr), 
                               region_to_load, refGen, verbose)   
        phased_snps_in_haploblock <-  as_tibble(loaded_vcf) %>%
          filter(str_detect(gt, "\\|")) %>%
          filter(!gt=="1|1") 
        if(nrow(phased_snps_in_haploblock)==0){
          #vm("no phased snps/heterozygous SNPS in haploblock", verbose)
          status_per_comb <- NULL
        } else {
          phased_snps_in_haploblock_ids <- phased_snps_in_haploblock %>%
            mutate(snp_id=paste0("s", c(1:nrow(.))))
          store_log(geneDir, phased_snps_in_haploblock_ids,
                    paste0("all_snps_hb_", 
                           HAP_ID, ".tsv"))
          ## get every mutation located in the haploblock
          muts_in_hap <- df_gene_hap %>%
            filter(hap_id==HAP_ID)
          ## calculate distance of each snp to each variant
          df_dist <- apply(muts_in_hap, 1, function(mut){
            phased_snps_in_haploblock_ids %>%
              mutate(dist=abs(start-as.numeric(mut[["pos"]])),
                     mut_id=mut[["mut_id"]]) %>%
              return()
          }) %>% bind_rows() %>%
            filter(dist<distCutOff) %>%
            arrange(dist) %>%
            mutate(comb_id=paste(mut_id, snp_id, sep="-"))
          store_log(geneDir, df_dist, paste0("dist_df_hb_", 
                                                         HAP_ID, ".tsv"))
          ## check if at least two main variants are present, otherwise phasing
          ## does not make any sense
          if(length(unique(df_dist$mut_id))>1){
            phased_per_mut_full <- phase_per_variant(df_dist, df_gene, verbose, 
                                                     bamDna, bamRna, purity,
                                                     printLog)
            phased_per_mut <- bind_rows(lapply(phased_per_mut_full, nth, 1))
            if(nrow(phased_per_mut)>0){
              log_phasing_results <- bind_rows(
                lapply(phased_per_mut_full, nth, 2)) %>%
                mutate(hap_id=HAP_ID)
              store_log(geneDir, log_phasing_results, paste0("log_phasing_hb_", 
                                                             HAP_ID, ".tsv"))
              store_log(geneDir, phased_per_mut, paste0("phased_per_mut_hb_", 
                                                        HAP_ID, ".tsv"))
              status_per_comb <- 
                define_status_for_main_combination(phased_per_mut, 
                                                   log_phasing_results)
            } else {
              status_per_comb <- NULL
            }
          } else {
            #vm("only one mut closer than dist limit to snp", verbose)
            status_per_comb <- NULL
          }
        }
        return(list(status_per_comb, log_phasing_results))
      })
      combined_status_per_haploblock <- lapply(per_haploblock, nth, 1) %>% 
        compact() %>%
        bind_rows() %>%
        select(comb, status, nstatus, conf) %>%
        mutate(phasing="haploblock")
      combined_phasing_per_haploblock <- lapply(per_haploblock, nth, 2) %>% 
        compact() %>%
        bind_rows() %>%
        mutate(phasing="haploblock")
      if(nrow(combined_status_per_haploblock)==0){
        combined_status_per_haploblock <- NULL
      } 
    } 
  } 
  return(list(status=combined_status_per_haploblock,
              info=combined_phasing_per_haploblock))
}