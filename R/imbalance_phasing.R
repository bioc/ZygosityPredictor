define_genotype_for_snps <- function(df_gene_seg, somCna, phasedVcf, distCutOff, 
                                     refGen, verbose){
  roi <- GRanges(paste0(unique(df_gene_seg$chr), ":",
                        min(df_gene_seg$pos), "-",
                        max(df_gene_seg$pos)))
  relevant_segments <- subsetByOverlaps(somCna, roi) 
  ## now load all snps which are inside the segment and inside the 
  ## distcutoff of the variants
  rel_snp_region <- df_gene_seg %>%
    mutate_at(.vars = c("pos"),
              .funs = as.numeric) %>%
    rowwise() %>%
    mutate(
      start=max(pos-distCutOff, start(relevant_segments)),
      end=min(pos+distCutOff, end(relevant_segments))
    ) %>%
    select(chr, start, end) %>%
    GRanges()
  loaded_snps <- loadVcf(phasedVcf, unique(df_gene_seg$chr), 
                         rel_snp_region, refGen, verbose)
  if(!"af" %in% nm_md(loaded_snps)){
    if("dp4" %in% nm_md(loaded_snps)){
      loaded_snps$af <- lapply(seq(1,length(loaded_snps)), function(i){
        sum(loaded_snps$dp4[[i]][c(3,4)])/sum(loaded_snps$dp4[[i]])
      }) %>% as.numeric()  
    } else {
      ## can not calculate AF
    }
  } 
  loaded_snps_hz <- loaded_snps[which(loaded_snps$af<0.95)]
  loaded_snps_hz$aff_cp <- lapply(loaded_snps_hz$af, aff_germ_copies, 
                                  chr=unique(df_gene_seg$chr), 
                                  tcn=unique(df_gene_seg$tcn),
                                  sex=sex,
                                  purity=purity) %>% unlist()
  loaded_snps_hz$rnd_aff_cp <- round(loaded_snps_hz$aff_cp)
  
  ## here the confidence measurement can be included
  loaded_snps_hz$gt <- case_when(
    loaded_snps_hz$rnd_aff_cp==min(loaded_snps_hz$rnd_aff_cp) ~ "1|0",
    loaded_snps_hz$rnd_aff_cp==max(loaded_snps_hz$rnd_aff_cp) ~ "0|1",
    TRUE ~ NA
  ) 
  snps <- as_tibble(loaded_snps_hz) %>% as_tibble() %>%
    mutate(snp_id=paste0("s", c(1:nrow(.))))
  return(snps)
}

perform_imbalance_phasing <- function(all_combinations, direct_phasing, 
                          haploblock_phasing, df_gene, somCna, 
                          phasedVcf, bamDna, bamRna, purity,  
                          distCutOff, printLog, 
                          verbose, geneDir, refGen){
  func_start(sys.call(), verbose)
  combined_phasing <- per_segment <- NULL
  do_imbalance_phasing <- decide_following_phasing(list(direct_phasing, 
                                                        haploblock_phasing), 
                                                   verbose)
  if(do_imbalance_phasing){
    df_gene_muts_all_imb <- df_gene %>% filter(all_imb==TRUE)
    if(nrow(df_gene_muts_all_imb)>1){
      df_gene_seg_full <- annotate_haploblocks_to_variants(df_gene_muts_all_imb, 
                                                           somCna, verbose)
      ## check if at least two variants are in the same haploblock
      if(any(df_gene_seg_full$seg_id > 0)&any(table(
        df_gene_seg_full[which(df_gene_seg_full$seg_id>0),]$seg_id
      )>1)){
        relevant_segments <- table(df_gene_seg_full$seg_id) %>% 
          .[which(.>1)] %>% names() %>% .[which(!.=="0")]
        per_segment <- lapply(relevant_segments, function(SEG){
          df_gene_seg <- df_gene_seg_full %>%
            filter(seg_id==SEG)
          snps <- define_genotype_for_snps(df_gene_seg, somCna, phasedVcf, 
                                           distCutOff, refGen, verbose)
          phasing_result <- phase_against_phased_snps(df_gene_seg, snps,  
                                                      verbose, bamDna, bamRna, 
                                                      purity, printLog, geneDir, 
                                                      distCutOff)
          return(phasing_result)
        })
        combined_phasing <- finalize_phasing_result(per_segment, "imbalance", 
                                                    all_combinations, verbose)
      } else {
        imbalance_phasing_exit <- "less than two variants in segment"
      }## not enough variants in same segment    
    } else {
      imbalance_phasing_exit <- "less than two variants have allelic imbalance"
    } ## not enough variants have alleleic imbalance
  }else {
    imbalance_phasing_exit <- "do imbalance = FALSE"
  } ## do_imbalance FALSE
  if(is.null(combined_phasing)){
    combined_phasing <- list(
      status=fill_phasing_status(all_combinations, NULL, 
                                 "imbalance", verbose),
      info=NULL,
      exit=imbalance_phasing_exit
    )
  }
  func_end(verbose)
  return(combined_phasing)
}


