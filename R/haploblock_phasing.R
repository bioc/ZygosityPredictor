#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
perform_haploblock_phasing <- function(df_gene, vcf, haploBlocks, phasedVcf,
                                    bamDna, bamRna, 
                                    purity, snp_dist=10000, distCutOff,
                                    printLog, verbose,
                                    geneDir, refGen){
  
  catt(printLog, 4, "Initializing level 2 phasing")
  #print(paste("GENE:", unique(df_gene$gene)))
  ## annotate to df_gene in which haploblock the variant is located, 0 means no 
  ## haploblock
  
  
  df_gene_hap <- apply(df_gene, 1, function(mut){
    gr_roi <- GenomicRanges::GRanges(paste0(mut[["chr"]], ":", 
                                            as.numeric(mut[["pos"]]), "-",
                                            as.numeric(mut[["pos"]])))
    #print(gr_roi)
    #print(haploBlocks)
    hap <- IRanges::subsetByOverlaps(haploBlocks, gr_roi)
    if(length(hap)==0){
      det_hap_id <- 0
    } else {
      det_hap_id <- hap$hap_id
      #print(mut)
      #print(hap)
    }
    return(as_tibble(t(mut)) %>% mutate(hap_id=det_hap_id))
  }) %>% bind_rows()
  
  #print(df_gene_hap)
  ## check if at least two variants are in the same haploblock
  if(any(df_gene_hap$hap_id > 0)&any(table(
    df_gene_hap[which(df_gene_hap$hap_id>0),]$hap_id
  )>1)){
    vm("at least two variants detected in same haploblock", verbose)
    #write_tsv(df_gene_hap, file=file.path("/home/m168r/home_extension/ZP_revision", paste0(unique(df_gene$gene), ".tsv")))
    relevant_haploblocks <- table(df_gene_hap$hap_id) %>% .[which(.>1)] %>% names() %>% .[which(!.=="0")]
    #print("relevant haploblocks")
    #print(relevant_haploblocks)
    ## iterate over haploblocks
    per_haploblock <- lapply(relevant_haploblocks, function(HAP_ID){
      #print(haploBlocks)
      #print(HAP_ID)
      region_to_load <- haploBlocks[which(haploBlocks$hap_id==HAP_ID)]
      
      # vcf_to_load <- phasedVcf %>% .[which(str_detect(., paste0("chr",unique(df_gene$chr), ".vcf")))]
      # #print(region_to_load)
      # #print("loading vcf")
      # #print(vcf_to_load)
      # vcf_region <- 
      #   VariantAnnotation::readVcf(vcf_to_load, "hg19")
      # 
      # #print("vcf loaded")
      # gt <- VariantAnnotation::geno(vcf_region)$GT %>% as.character()
      # comb_vcf <- rowRanges(vcf_region)
      # comb_vcf$gt <- gt
      # 
      # #print(comb_vcf)
      # #print(region_to_load)
      # 
      loaded_vcf <-  loadVcf(phasedVcf, unique(df_gene$chr), 
                             region_to_load, refGen, verbose) #%>%
      #vm(loaded_vcf, verbose)   
      phased_snps_in_haploblock <-  as_tibble(loaded_vcf) %>%
        filter(str_detect(gt, "\\|")) %>%
        filter(!gt=="1|1") 
      
      if(nrow(phased_snps_in_haploblock)==0){
        vm("no phased snps/heterozygous SNPS in haploblock", verbose)
        return(NULL)
      } else {
        phased_snps_in_haploblock_ids <- phased_snps_in_haploblock %>%
          mutate(snp_id=paste0("s", c(1:nrow(.))))
        if(!is.null(geneDir)){
          write_tsv(phased_snps_in_haploblock_ids, 
                    file=file.path(geneDir, 
                                   paste0("all_snps_hb_", 
                                          HAP_ID, ".tsv")))
        }
        #print(phased_snps_in_haploblock)
        
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
        
        #print(df_dist)
        if(!is.null(geneDir)){
          write_tsv(df_dist, 
                    file=file.path(geneDir, 
                                   paste0("full_dist_df_hb_", 
                                          HAP_ID, ".tsv")))
        }
        if(length(unique(df_dist$mut_id))>1){
          phased_per_mut_full <- lapply(unique(df_dist$mut_id), function(MUT_ID){
            #message("phasing", MUT_ID)
            df_dist_mut <- df_dist %>%
              filter(mut_id==MUT_ID)
            #print(df_dist_mut)
            high_conf_result <- FALSE
            i <- 1
            phased_muts <- tibble()
            ## only for documentation
            sub_status_list <- tibble()
            while(high_conf_result==FALSE&i<=nrow(df_dist_mut)){
              cat(i)
              
              comb <- df_dist_mut[i,]
              #print("hier kommt comb")
              #print(comb)
              sub_comb <- bind_rows(
                df_gene %>% 
                  filter(mut_id==comb$mut_id) %>%
                  select(chr, pos, mut_id, ref, alt, class),
                comb %>%
                  mutate(ALT=unlist(
                    lapply(ALT, function(x){as.character(x[[1]][1])})),
                    class="snv") %>%
                  select(chr=seqnames, pos=start, mut_id=snp_id, ref=REF, alt=ALT, 
                         class)
              ) %>%
                make_phasing_combinations()
              
              ## the function classify_reads is usually used inside of apply()
              ## thats why the input needs to be adjusted in the format of apply iterations:
              comb_vec <- as.character(sub_comb) %>% set_names(nm=names(sub_comb))
              
              
              
              res <- classify_reads(comb_vec, bamDna, bamRna, verbose)
              
              if(nrow(res)!=0){
                sub_status_combination <- 
                  classify_combination(res,
                                       purity,
                                       eval_full=FALSE,
                                       printLog,
                                       verbose)
                #vm("111111111111", verbose)
                sub_status_list <- bind_rows(sub_status_list,
                                             sub_status_combination %>% 
                                               mutate(comb_id=comb$comb_id))
                #print(sub_status_list)
                #vm("222222221", verbose)
                if(!sub_status_combination$status=="null"){
                  phased <- tibble(mut_id=comb$mut_id,
                                   gt=get_genotype(comb$gt, sub_status_combination$status))
                  phased_muts <- bind_rows(phased_muts, phased)
                  
                  high_conf_result <- TRUE
                  
                } 
                
              }
              i <- i +1 
              
            } 
            
            
            return(list(phased_muts, sub_status_list))
            
            
            # if(length(bam_raw)>0){
            #   write_tsv(tibble(x="hello"), file=file.path("/home/m168r/home_extension/ZP_revision", paste0(unique(df_gene$gene), ".tsv")))
            # }
            
          }) #%>% bind_rows()
          #vm("333333333", verbose)
          phased_per_mut <- bind_rows(lapply(phased_per_mut_full, nth, 1))
          
          #print("done phasing of variants to snps")
          #print(phased_per_mut)
          if(nrow(phased_per_mut)>0){
            #print(phased_per_mut_full)
            if(!is.null(geneDir)){
              log_phasing_results <- bind_rows(lapply(phased_per_mut_full, nth, 2))
              write_tsv(log_phasing_results, 
                        file=file.path(geneDir, 
                                       paste0("phasing_info_hb_", 
                                              HAP_ID, ".tsv")))
            }
            
            #vm("111111111111", verbose)
            status_per_comb <- outer(phased_per_mut$mut_id, 
                                     phased_per_mut$mut_id, `paste`) %>%
              .[which(upper.tri(.))] %>%
              as.data.frame() %>% 
              set_names(nm='raw') %>%
              mutate(
                mut_id1=str_split(raw, " ") %>% map_chr(.,1),
                mut_id2=str_split(raw, " ") %>% map_chr(.,2)
              ) %>%
              left_join(phased_per_mut %>% select(mut_id1=mut_id, gt1=gt), 
                        by="mut_id1") %>%
              left_join(phased_per_mut %>% select(mut_id2=mut_id, gt2=gt), 
                        by="mut_id2") %>%
              rowwise() %>%
              mutate(status=case_when(gt1==gt2 ~ "same",
                                      TRUE ~ "diff"),
                     #comb_id=paste(mut_id1, mut_id2, sep='-')
                     comb_id=paste(paste0("m",
                                          sort(as.numeric(str_match(
                                            c(mut_id1, mut_id2), "\\d+")))),
                                   collapse = "-")
              ) %>%
              ungroup() %>%
              select(comb_id, status)  
          } else {
            status_per_comb <- NULL
          }
          
          #vm("222222222", verbose)
          return(status_per_comb)
          #print(table(phased_per_mut))
          
          
          # if(nrow(phased_per_mut)>1){
          #   write_tsv(tibble(x="hello"), file=file.path("/home/m168r/home_extension/ZP_revision/gene_out", paste0(unique(df_gene$gene), ".tsv")))
          # }        
        } else {
          vm("only one mut closer than dist limit to snp", verbose)
          return(NULL)
        }
      }
    })  %>%
      bind_rows()
    
    if(nrow(per_haploblock)==0){
      return(NULL)
    } else {
      return(per_haploblock)
    }
  } else {
    #print("no pair of variants located in the same haploblock")
    return(NULL)
  }
}