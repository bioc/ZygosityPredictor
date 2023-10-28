define_genotype_for_snps <- function(df_gene_seg, somCna, phasedVcf, distCutOff, 
                                     refGen, purity, sex, verbose){
  func_start(verbose)
  ##prind(1)
  roi <- GRanges(paste0(unique(df_gene_seg$chr), ":",
                        min(df_gene_seg$pos), "-",
                        max(df_gene_seg$pos)))
  ##prind(roi)
  relevant_segments <- subsetByOverlaps(somCna, roi) 
  major_allele <- max(as.numeric(unlist(str_split(relevant_segments$genotype, ":"))))
  minor_allele <- min(as.numeric(unlist(str_split(relevant_segments$genotype, ":"))))
  ## now load all snps which are inside the segment and inside the 
  ## distcutoff of the variants
  ##prind(2)
  # #prind(relevant_segments)
  # #prind(end(relevant_segments))
  # #prind(GenomicRanges::end(relevant_segments))
  # #prind(start(relevant_segments)) 
  df_gene_seg_num <- df_gene_seg %>%
    mutate_at(.vars = c("pos"),
              .funs = as.numeric) %>%
    rowwise() #%>%
  reg_upstream <- df_gene_seg_num %>%
    mutate(
      ## there is another start and end function, --> GenomicRanges::
      start=max(pos-distCutOff, GenomicRanges::start(relevant_segments)),
      end=pos-1
    ) %>%
    select(chr, start, end)
  reg_downstream <- df_gene_seg_num %>%
    mutate(
      ## there is another start and end function, --> GenomicRanges::
      start=pos+1,
      end=min(pos+distCutOff, GenomicRanges::end(relevant_segments))
    ) %>%
    select(chr, start, end)
  rel_snp_region_pre <- bind_rows(reg_upstream, reg_downstream)  
  #prind(rel_snp_region_pre)
  rel_snp_region <- GRanges(rel_snp_region_pre)
  ##prind(0, debug)
  ##prind(phasedVcf, debug)
  ##prind(3)
  ##prind(phasedVcf)
  loaded_snps <- loadVcf(phasedVcf, unique(df_gene_seg$chr), 
                         rel_snp_region, refGen, verbose)
  ##prind(4)
  ##prind(loaded_snps)
  ##prind(loaded_snps$gt)
  if(!is.null(loaded_snps)){
    if(!"af" %in% nm_md(loaded_snps)){
      if("dp4" %in% nm_md(loaded_snps)){
        loaded_snps$af <- lapply(seq(1,length(loaded_snps)), function(i){
          sum(loaded_snps$dp4[[i]][c(3,4)])/sum(loaded_snps$dp4[[i]])
        }) %>% as.numeric()  
      } 
    } 
    if("af" %in% nm_md(loaded_snps)){
      ##prind(5)
      ##prind(elementMetadata(loaded_snps))
      loaded_snps_hz <- loaded_snps[which(loaded_snps$af<0.95&str_detect(loaded_snps$gt, "0"))]
      if(length(loaded_snps_hz)>0){
        ##prind(loaded_snps_hz$af)
        ##prind(loaded_snps_hz$dp4)
        ##prind(loaded_snps_hz["dp4"])
        loaded_snps_hz$aff_cp <- lapply(loaded_snps_hz$af, aff_germ_copies, 
                                        chr=unique(df_gene_seg$chr), 
                                        tcn=unique(df_gene_seg$tcn),
                                        sex=sex,
                                        purity=purity) %>% unlist()
        ##prind(loaded_snps_hz$aff_cp)
        loaded_snps_hz$rnd_aff_cp <- round(loaded_snps_hz$aff_cp)
        ##prind(3, debug)
        ## here the confidence measurement can be included
        ##prind(6)
        ##prind(loaded_snps_hz)
        if(length(loaded_snps_hz)==1){
          ## if there is only one it doesnt matter
          loaded_snps_hz$gt <- "0|1"
        } else {
          otsu <- calc_otsu(loaded_snps_hz$aff_cp)
          loaded_snps_hz$gt <- case_when(
            loaded_snps_hz$aff_cp < otsu ~ "1|0",
            loaded_snps_hz$aff_cp > otsu ~ "0|1",
            #abs(loaded_snps_hz$rnd_aff_cp-minor_allele)<abs(loaded_snps_hz$rnd_aff_cp-major_allele) ~ "1|0",
            #abs(loaded_snps_hz$rnd_aff_cp-minor_allele)>abs(loaded_snps_hz$rnd_aff_cp-major_allele) ~ "0|1",
            # loaded_snps_hz$rnd_aff_cp==min(loaded_snps_hz$rnd_aff_cp) ~ "1|0",
            # loaded_snps_hz$rnd_aff_cp==max(loaded_snps_hz$rnd_aff_cp) ~ "0|1",
            TRUE ~ NA
          )         
        }
  
        ##prind(7)
        snps <- as_tibble(loaded_snps_hz) %>%
          mutate(snp_id=paste0("s", c(1:nrow(.))))
        append_loglist(nrow(snps),"heterozygous SNPs detected")
        ##prind("snps")
        ##prind(snps)        
      } else {
        ## only homozygous snps
        snps <- NULL
        append_loglist("Only homozygous SNPs")
      }

    } else {
      ## AF can not be calculated
      snps <- NULL
      append_loglist("no annoatted Allele-frequency/ no calculation from dp4 possible")
    }
    ##prind(2, debug)

  } else {
    snps <- NULL
    append_loglist("no vcf file or no SNPS in roi")
    ## no vcf / no snps in roi
  }
  func_end(verbose)
  return(snps)
}

perform_imbalance_phasing <- function(all_combinations, direct_phasing, 
                          haploblock_phasing, df_gene, somCna, 
                          phasedVcf, bamDna, bamRna, purity, sex,  
                          distCutOff, printLog, 
                          verbose, geneDir, refGen){
  func_start(verbose)
  phasing_type <- "imbalance"
  combined_phasing <- per_segment <- phasingDir <- NULL
  ## pre define empty result
  combined_phasing <- list(
    status=fill_phasing_status(all_combinations, NULL, 
                               phasing_type, verbose),
    info=NULL)
  do_imbalance_phasing <- decide_following_phasing(list(direct_phasing, 
                                                        haploblock_phasing), 
                                                   phasedVcf,
                                                   verbose)
  if(do_imbalance_phasing){
    append_loglist("Initializing Imbalance phasing")
    if(!is.null(geneDir)){
      phasingDir <- file.path(geneDir, phasing_type)
      dir.create(phasingDir)
    }
    df_gene_muts_all_imb <- df_gene %>% filter(all_imb==TRUE)
    if(nrow(df_gene_muts_all_imb)>1){
      df_gene_seg_full <- annotate_haploblocks_to_variants(df_gene_muts_all_imb, 
                                                           somCna, verbose)
      ## check if at least two variants are in the same haploblock
      if(any(df_gene_seg_full$seg_id > 0)&any(table(
        df_gene_seg_full[which(df_gene_seg_full$seg_id>0),]$seg_id
      )>1)){
        append_loglist("At least 2 variants in segment of allelic imbalance")
        relevant_segments <- table(df_gene_seg_full$seg_id) %>% 
          .[which(.>1)] %>% names() %>% .[which(!.=="0")]
        per_segment <- lapply(relevant_segments, function(SEG){
          df_gene_seg <- df_gene_seg_full %>%
            filter(seg_id==SEG)
          snps <- define_genotype_for_snps(df_gene_seg, somCna, phasedVcf, 
                                           distCutOff, refGen, purity, sex, 
                                           verbose)
          store_log(phasingDir, snps, "snps_check_genotype_assignment.tsv")
          if(!is.null(snps)){
            phasing_result <- phase_against_phased_snps(df_gene_seg, snps,  
                                                        verbose, bamDna, bamRna, 
                                                        purity, printLog, 
                                                        phasingDir, 
                                                        phasing_type,
                                                        distCutOff)            
          } else {
            phasing_result <- list(NULL, NULL, 
                                   "no snps/ no vcf in/for roi")
          }

          return(phasing_result)
        })
        combined_phasing <- finalize_phasing_result(per_segment, phasing_type, 
                                                    all_combinations, verbose)
        imbalance_phasing_exit <- combined_phasing[[3]]
      } else {
        imbalance_phasing_exit <- "less than two variants in segment"
      }
    } else {
      imbalance_phasing_exit <- "less than two variants have allelic imbalance"
    } 
  } else {
    imbalance_phasing_exit <- "do imbalance = FALSE"
  } 
  append_loglist(imbalance_phasing_exit)
  if(nrow(combined_phasing$status)==0){
    combined_phasing$status <- fill_phasing_status(all_combinations, NULL, 
                                                   phasing_type, verbose)
  }
  full_phasing <- append(combined_phasing[c(1,2)], 
                         list(exit=imbalance_phasing_exit))
  
  func_end(verbose)
  #print(combined_phasing)
  return(full_phasing)
}


