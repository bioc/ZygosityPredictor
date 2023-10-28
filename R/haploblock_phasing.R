phase_per_variant <- function(df_dist, muts_in_hap, verbose, bamDna, bamRna, 
                              purity,
                              printLog, phasing_type){
  func_start(verbose)
  ##prind(df_dist)
  #prind(0)
  phased_per_mut_full <- lapply(unique(df_dist$mut_id), function(MUT_ID){
    
    df_dist_mut <- df_dist %>%
                     filter(mut_id==MUT_ID)
    ##prind(df_dist_mut)
    high_conf_result <- FALSE
    i <- 1
    phased_muts <- tibble()
    ## only for documentation
    sub_status_list <- tibble()
    while(high_conf_result==FALSE&i<=nrow(df_dist_mut)){
      #prind(1)
      comb <- df_dist_mut[i,]
      #prind(MUT_ID)
      #prind(comb)
      #prind(muts_in_hap)
      sub_comb <- bind_rows(
                    muts_in_hap %>% 
                      filter(mut_id==comb$mut_id) %>%
                      select(chr, pos, mut_id, ref, alt, class) %>%
                      mutate_all(.funs = as.character),
                    comb %>%
                      mutate(ALT=unlist(
                        lapply(ALT, function(x){as.character(x[[1]][1])})),
                          class="snv") %>%
                          select(chr=seqnames, pos=start, mut_id=snp_id, 
                                 ref=REF, alt=ALT, class) %>%
                      mutate_all(.funs = as.character)
                        ) %>%
        make_phasing_combinations() %>%
        ## THIS IS VERY IMPORTANT; FOR SOME REASON THE TRANSITION FROM FACTOR TO CHARACTERE IN THE NEXT FUNCTION CAUSES A -1 FOR THE CHROMOSOMSE!!!!!!
        mutate(chr1=as.character(chr1), chr2=as.character(chr2))
      append_loglist(phasing_type, "phasing of", sub_comb[["comb_id"]], "distance:", sub_comb[["dist"]])
      #prind(sub_comb)
 ## the function classify_reads is usually used inside of apply()
## thats why the input needs to be adjusted in the format of apply iterations:
      comb_vec <- as.character(sub_comb) %>% set_names(nm=names(sub_comb))
      #prind(comb_vec)
      
      res <- classify_reads(comb_vec, bamDna, bamRna, verbose)
      #prind(2)
      if(nrow(res)!=0){
        #prind(3)
        #prind(res)
        append_loglist(nrow(res), 
                       "reads / read-pairs covering both positions")
        sub_status_combination <- classify_combination(res,
                                                       purity,
                                                       printLog,
                                                       verbose) %>%
          mutate(dist=comb$dist,
                 class_comb=
                   paste(muts_in_hap[which(muts_in_hap$mut_id==comb$mut_id),]$class,
                         "snp",
                         sep="-"))
        sub_status_list <- bind_rows(sub_status_list,
                                     sub_status_combination %>% 
                                       mutate(comb=comb$comb_id))
        ##prind(sub_status_list)
        #prind(1)
        if(!sub_status_combination$status=="null"){
          #prind(comb)
          
          phased <- tibble(mut_id=comb$mut_id,
                           gt=get_genotype(comb$gt, 
                             sub_status_combination$status))
          phased_muts <- bind_rows(phased_muts, phased)
          high_conf_result <- TRUE
          local_haploblock_phasing_exit <- paste0(comb$mut_id, 
                                                  ": completed")
        } else {
          local_haploblock_phasing_exit <- paste0(comb$mut_id, 
                                                  ": phasing status null")
        }
      } else {
        local_haploblock_phasing_exit <- paste0(comb$mut_id, 
                                                ": no covering reads")
      }
      i <- i +1 
      append_loglist(local_haploblock_phasing_exit)
    } 
    return(list(phased_muts, sub_status_list, local_haploblock_phasing_exit)) 
  }) 
  func_end(verbose)
  return(phased_per_mut_full)
}
annotate_haploblocks_to_variants <- function(df_gene, gr_obj, verbose){
  func_start(verbose)
  ## annotate to df_gene in which haploblock the variant is located, 0 means no 
  ## haploblock
  ## first check if object is somCna or haploBlocks
  ##prind(1,debug)
  id_col <- str_match(nm_md(gr_obj), "hap_id|seg_id") %>% .[which(!is.na(.))]
  ##prind(1,debug)
  df_gene_hap <- apply(df_gene, 1, function(mut){
    gr_roi <- GenomicRanges::GRanges(paste0(mut[["chr"]], ":", 
                                            as.numeric(mut[["pos"]]), "-",
                                            as.numeric(mut[["pos"]])))
    hap <- IRanges::subsetByOverlaps(gr_obj, gr_roi)
    if(length(hap)==0){
      det_hap_id <- 0
    } else {
      det_hap_id <- elementMetadata(hap)[[id_col]]
    }
    return(as_tibble(t(mut)) %>% mutate(!!id_col:=det_hap_id))
  }) %>% bind_rows()
  ##prind(1,debug)
  func_end(verbose)
  return(df_gene_hap)
}
define_status_for_main_combination <- function(phased_per_mut, 
                                               log_phasing_results,
                                               verbose){
  vm(as.character(sys.call()[1]), verbose, 1)
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
  vm("  - done", verbose, -1)
  return(status_per_comb)
}

phase_against_phased_snps <- function(muts_in_hap, 
                                      phased_snps_in_haploblock_ids, 
                                      verbose, 
                                      bamDna, bamRna, purity,
                                      printLog, phasingDir, phasing_type, distCutOff){
  func_start(verbose)
  name_id_col <- str_match(names(muts_in_hap), "hap_id|seg_id") %>% .[!is.na(.)]
  HAP_ID <- unique(muts_in_hap[[name_id_col]])
  log_phasing_results <- NULL
  df_dist <- apply(muts_in_hap, 1, function(mut){
    phased_snps_in_haploblock_ids %>%
      mutate(dist=abs(start-as.numeric(mut[["pos"]])),
             mut_id=mut[["mut_id"]]) %>%
      return()
  }) %>% bind_rows() %>%
    filter(dist<distCutOff) %>%
    arrange(dist) %>%
    mutate(comb_id=paste(mut_id, snp_id, sep="-"))
 ##prind(df_dist)
  store_log(phasingDir, df_dist, paste0("dist_df_hb_", 
                                     HAP_ID, ".tsv"))
  ## check if at least two main variants are present, otherwise phasing
  ## does not make any sense
  if(length(unique(df_dist$mut_id))>1){
    append_loglist("at least 2 variants closer to SNPs than distCutOff")
    phased_per_mut_full <- phase_per_variant(df_dist, muts_in_hap, verbose, 
                                             bamDna, bamRna, purity,
                                             printLog, phasing_type)
    ##prind(phased_per_mut_full)
    phased_per_mut <- bind_rows(lapply(phased_per_mut_full, nth, 1))
    local_haploblock_phasing_exit <- lapply(phased_per_mut_full, nth, 3) %>% 
      unlist() %>% paste(collapse = "; ")
    if(nrow(phased_per_mut)>0){
      ## log_phasing_result is a tibble of all combinations with main muts and
      ## snps in the haploblock... can be df with many rows
      log_phasing_results <- bind_rows(
        lapply(phased_per_mut_full, nth, 2)) %>%
        mutate(hap_id=HAP_ID)
      
      store_log(phasingDir, log_phasing_results, paste0("log_phasing_hb_", 
                                                     HAP_ID, ".tsv"))
      store_log(phasingDir, phased_per_mut, paste0("phased_per_mut_hb_", 
                                                HAP_ID, ".tsv"))
      status_per_comb <- 
        define_status_for_main_combination(phased_per_mut, 
                                           log_phasing_results, verbose)
    } else {
      status_per_comb <- NULL
    }
  } else {
    local_haploblock_phasing_exit <- 
      "only one mut closer than distCutOff to SNP"
    append_loglist(local_haploblock_phasing_exit)
    #vm("only one mut closer than dist limit to snp", verbose)
    status_per_comb <- NULL
  }
  func_end(verbose)
  return(list(status_per_comb, log_phasing_results, 
              local_haploblock_phasing_exit))
}

#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
perform_haploblock_phasing <- function(all_combinations, direct_phasing, df_gene, 
                                       haploBlocks, phasedVcf,
                                    bamDna, bamRna, 
                                    purity, distCutOff,
                                    printLog, verbose,
                                    geneDir, refGen){
  func_start(verbose)
  per_haploblock <- log_phasing_results <- NULL
  phasing_type <- "haploblock"
  ## pre define empty result
  combined_phasing <- list(
    status=fill_phasing_status(all_combinations, NULL, 
                               phasing_type, verbose),
    info=NULL)
  
  do_haploblock_phasing <- eval_direct_results(direct_phasing, phasedVcf,
                                               haploBlocks)
  if(do_haploblock_phasing){
    append_loglist("Initializing haploblock phasing")
    if(!is.null(geneDir)){
      phasingDir <- file.path(geneDir, phasing_type)
      dir.create(phasingDir)
    }
    df_gene_hap <- annotate_haploblocks_to_variants(df_gene, haploBlocks, 
                                                    verbose)
    ## check if at least two variants are in the same haploblock
    if(any(df_gene_hap$hap_id > 0)&any(table(
      df_gene_hap[which(df_gene_hap$hap_id>0),]$hap_id
    )>1)){
      append_loglist("at least 2 variants detected in same haploblock")
      #vm("at least two variants detected in same haploblock", verbose, 1)
      relevant_haploblocks <- table(df_gene_hap$hap_id) %>% 
        .[which(.>1)] %>% names() %>% .[which(!.=="0")]
      ## iterate over haploblocks
      per_haploblock <- lapply(relevant_haploblocks, function(HAP_ID){
        append_loglist("checking haploblock:", HAP_ID)
        log_phasing_results <- NULL
        region_to_load <- haploBlocks[which(haploBlocks$hap_id==HAP_ID)]
        loaded_vcf <-  loadVcf(phasedVcf, unique(df_gene$chr), 
                               region_to_load, refGen, verbose) 
        ###prind(loaded_vcf, debug)
        append_loglist(length(loaded_vcf), "SNPs in haploblock")
        phased_snps_in_haploblock <-  as_tibble(loaded_vcf) %>%
          ## 
          filter(str_detect(gt, "\\|")) %>%
          filter(!gt=="1|1") %>%
          filter(!as.numeric(start) %in% df_gene$pos)
        append_loglist(nrow(phased_snps_in_haploblock), "phased &  heterozygous")
        #prind(phased_snps_in_haploblock)
        if(!nrow(phased_snps_in_haploblock)==0){
          phased_snps_in_haploblock_ids <- phased_snps_in_haploblock %>%
            mutate(snp_id=paste0("s", c(1:nrow(.))))
          #print(phased_snps_in_haploblock_ids)
          store_log(phasingDir, phased_snps_in_haploblock_ids,
                    paste0("all_snps_hb_", 
                           HAP_ID, ".tsv"))
          ## get every mutation located in the haploblock
          muts_in_hap <- df_gene_hap %>%
            filter(hap_id==HAP_ID)
          phasing_result <- phase_against_phased_snps(muts_in_hap, 
                                    phased_snps_in_haploblock_ids, 
                                                verbose, 
                                                bamDna, bamRna, purity,
                                                printLog, phasingDir, 
                                    phasing_type, distCutOff)         #vm("no phased snps/heterozygous SNPS in haploblock", verbose)
          ###prind(phasing_result, debug)
          #local_haploblock_phasing_exit <- phasing_result[[3]]
        } else {
          ## hier empty
          local_haploblock_phasing_exit <- 
            paste("no phased snps in haploblock:", HAP_ID)
          phasing_result <- list(tibble(), NULL, local_haploblock_phasing_exit)
        }
        return(phasing_result)
      })
      ##prind(per_haploblock, debug)
      combined_phasing <- finalize_phasing_result(per_haploblock, phasing_type, 
                                                  all_combinations, verbose)
      global_haploblock_phasing_exit <- combined_phasing[[3]]
    } else {
      global_haploblock_phasing_exit <- "less than two variants in same haploblock"
    }
  } else {
    global_haploblock_phasing_exit <- "do_hap_phasing = FALSE"
  }
  append_loglist(global_haploblock_phasing_exit)
  #print(combined_phasing$status)
  if(nrow(combined_phasing$status)==0){
    combined_phasing$status <- fill_phasing_status(all_combinations, NULL, 
                                 phasing_type, verbose)
  }
  full_phasing <- append(combined_phasing[c(1,2)], 
                         list(exit=global_haploblock_phasing_exit))
  func_end(verbose)
  return(full_phasing)
}

finalize_phasing_result <- function(per_hapseg, phasing_type, all_combinations,
                                    verbose){
  func_start(verbose)
  combined_phasing_per_hapseg <- NULL
  combined_status_per_hapseg <- fill_phasing_status(all_combinations, 
                                                      NULL, 
                                                      phasing_type, verbose)
  if(length(compact(per_hapseg))>0){  
    ##prind(2,debug)  
    status_per_hapseg <- lapply(per_hapseg, nth, 1) %>% 
      compact() 
    if(length(status_per_hapseg)>0){
      combined_status_per_hapseg <- status_per_hapseg %>%
        bind_rows() %>%
        select(comb, status, nstatus, conf) %>%
        mutate(phasing=phasing_type) 
   
    } 
    phasing_per_hapseg <- lapply(per_hapseg, nth, 2) %>% 
        compact()       
    if(length(phasing_per_hapseg)>0){ 
      ###prind(5,debug)
      combined_phasing_per_hapseg <- phasing_per_hapseg %>%
        bind_rows() %>%
        mutate(phasing=phasing_type)    
    }  
  } 
  ##prind(4,debug)

  ###prind(phasing_per_hapseg, debug)
  exit_per_hapseg <- lapply(per_hapseg, nth, 3) %>% 
    unlist() %>% paste(collapse = "; ")
  
  ##prind(6,debug)
  func_end(verbose)
  return(list(status=combined_status_per_hapseg,
            info=combined_phasing_per_hapseg,
            exit=tibble(type=phasing_type,
                        info=exit_per_hapseg)))
}


