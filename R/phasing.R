get_single_index <- function(x,y,nr){
  nr*(x-1)+y
}
get_xy_index <- function(inp, nr){
  lapply(inp, function(i){
    y <- i%%nr
    if(y==0){
      y_ret <- nr
      x <- (i-y)/nr 
    } else {
      y_ret <- y
      x <- (i-y)/nr+1
    }
    return(c(x=x, y=y_ret))  
  }) %>% Reduce(function(x,y)rbind(x,y),.)
  
}
get_main_mut_pos <- function(){
  lm <- length(main_muts)
  lapply(1:lm, function(i){
    s <- (i-1)*nrow(mat_phased)+1
    e <- (i-1)*nrow(mat_phased)+lm
    return(seq(s,e,1))
  }) %>% unlist() %>%
    return()
}
unphased_main <- function(){
  lm <- seq(1,length(main_muts),1)
  which(is.na(mat_phased[lm, lm]))
}
unknown_main <- function(){
  all_unknown <- which(upper.tri(mat_phased)&mat_phased==0)
  intersect(all_unknown, get_main_mut_pos())
}
#' @importFrom stringr %>%
#' @importFrom dplyr select filter mutate bind_rows
#' @importFrom tibble as_tibble
add_snps_to_matrices <- function(snps){
  func_start()
  snps_hap <- snps[which(!is.na(snps$block_final)),] %>%
    select(mut_id, block_final, gt=gt_final)
  known_combs <- lapply(unique(snps_hap$block_final) %>% .[which(.!=0)], function(HAP_ID){
    
    snps_in_hap <- snps_hap %>%
      filter(block_final==HAP_ID) %>%
      mutate(id=as.numeric(str_match(mut_id, "\\d+")))
    if(nrow(snps_in_hap)>1){
      ## at least two snps in haploblock required to draw any concludios
      res_list <- list()
      for (i in 1:(nrow(snps_in_hap)-1)) {
        result_matrix <- matrix(0, nrow = nrow(snps_in_hap)-i, ncol = 2)
        for (j in (i+1):nrow(snps_in_hap)) {
          comparison_result <- ifelse(snps_in_hap$gt[i] == snps_in_hap$gt[j], 1, 2)
          # Store the results in the result_matrix
          result_matrix[nrow(snps_in_hap)-j+1, 1] <- get_single_index(snps_in_hap$id[j], snps_in_hap$id[i], nrow(mat_phased))
          result_matrix[nrow(snps_in_hap)-j+1, 2] <- comparison_result
        }
        res_list[[i]] <- result_matrix
      }
      to_add <- Reduce(function(x,y)rbind(x,y),res_list)
      colnames(to_add) <- c("comb", "nstatus") 
      new <- to_add %>%
        as_tibble() %>%
        mutate(conf=5, phasing=HAP_ID, hap_id=HAP_ID)     
    } else {
      ## if haploblock has only one snp, no concludions can be drawn
      new <- tibble()
    }
    return(new)
  })  %>% bind_rows()
  append_matrices(known_combs, iterate=FALSE)
  func_end()
}
eval_rare_case <- function(all_comb){
  diff_combs <- all_comb[which(all_comb$nstatus==2),]
  if(nrow(diff_combs)>2){
    ## more than 2 combinatiosna re diff
    ## check for all posiible master combinations
    mut_occurence <- unlist(str_split(diff_combs$comb, "-")) %>% table() %>%
      .[which(.>=2)]
    if(length(mut_occurence)>=3){
      rare_case <- " rare_case_detected"
      
      
      
    } else {
      rare_case <- ""
    }
  } else {
    rare_case <- ""
  }
  return(rare_case)
}
create_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  mat_dist <<- make_dist_matrix(all_pos, all_variants, distCutOff) 
  main_muts <<- all_variants
  main_pos <<- all_pos
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA  
  mat_info <<- mat_phased
}
append_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  func_start()
  mat_dist <<- make_dist_matrix(c(main_pos, all_pos), 
                                c(main_muts, all_variants), 
                                distCutOff) 
  mat_phased_main <- mat_phased
  mat_info_main <- mat_info
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA 
  mat_info <<- mat_phased
  nm <- length(main_muts)
  mat_phased[c(1:nm), c(1:nm)] <<- mat_phased_main
  mat_info[c(1:nm), c(1:nm)] <<- mat_info_main
  func_end()
}
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
get_main_mut_conns <- function(){
  lapply(main_muts, function(M){
    length(c(which(mat_phased[M,]!=0|is.na(mat_phased[M,])),
             which(mat_phased[,M]!=0|is.na(mat_phased[,M]))
    ))
  }) %>% unlist() %>%
    set_names(main_muts)
}
#' @importFrom magrittr %>%
#' @importFrom tibble as_tibble
#' @importFrom purrr set_names
#' @importFrom igraph graph_from_data_frame all_shortest_paths
get_next_path <- function(comb, distCutOff){
  func_start()
      
  mains <- as.character(sort(get_xy_index(as.numeric(comb), nrow(mat_phased))))

  ldc <- 1000
        open_conns <- which((is.na(mat_phased)|mat_phased>0)&mat_dist<ldc) %>%
          get_xy_index(.,nrow(mat_phased))
        if(length(dim(open_conns))!=2){
          connections <- as_tibble(t(open_conns)) %>%
            set_names(c("mut_id1", "mut_id2"))
        } else {
          colnames(open_conns) <- c("mut_id1", "mut_id2")
          connections <- as_tibble(open_conns)
        }
  if(length(intersect(as.numeric(mains), 
                      c(connections$mut_id1, connections$mut_id2)))==2){
    graph <- graph_from_data_frame(connections, directed = FALSE)
    all_paths <- igraph::all_shortest_paths(graph, from = mains[1], to = mains[2], mode = "all")$res  
    if(length(all_paths)>0){
        poss <- list()
        for (PATH in all_paths){
          P <- as.numeric(names(PATH))
          res <- c()
          for(i in c(1:(length(P)-1))){
            mv <- sort(c(P[i], P[i+1]))
            conn <- get_single_index(mv[2], mv[1], nrow(mat_phased))
            res[i] <- c(res, mat_phased[conn])
            if(is.na(mat_phased[conn])){
              poss <- c(poss, conn)
            }
          }
        } 
        np <-unlist(poss)[1]
      } else {
        np <- NULL
      } 
  } else {
    np <- NULL
  }
  func_end()
  return(np)
}
#' @importFrom purrr set_names
#' @importFrom stringr %>%
#' @importFrom dplyr left_join
prioritize_combination <- function(){
  func_start()
  unrel_rows <- rownames(mat_phased) %>% .[which(!. %in% colnames(mat_phased))]
  mat_tmp <- mat_dist
  mat_tmp[unrel_rows,] <- NA
  mat_tmp[,unrel_rows] <- NA
  mat_tmp[c((length(main_muts)+1):nrow(mat_phased)),] <- NA
  unphased <- which(is.na(mat_phased))
  relevant_dists <- mat_tmp[unphased] 
  relevant_dists[which(relevant_dists==0)] <- NA
  ## to suppress warning: Warning: no non-missing arguments to min; returning Inf
  suppressWarnings(
    shortest <- unphased[which(relevant_dists==min(relevant_dists, na.rm=T))] 
  )
  ## pick first one, if two have the same distance
  if(length(shortest)==0){
    main_mut_conns <- get_main_mut_conns()
    df_mmc <- tibble(m=names(main_mut_conns),
                     n=main_mut_conns)
    uc <- unknown_main()
    if(length(uc)>0){
      df_unknown_pre <- lapply(uc, function(C){
        sort(get_xy_index(C, nrow(mat_phased))) %>%
          paste0("m",.) %>%
          set_names(c("m1", "m2")) %>%
          c(.,comb=C)
      }) %>%
        bind_rows() 
      df_unknown <- df_unknown_pre %>%
        left_join(df_mmc %>% select(m, n1=n), by=c("m1"="m")) %>%
        left_join(df_mmc %>% select(m, n2=n), by=c("m2"="m")) %>%
        filter(n1>0&n2>0)
      if(nrow(df_unknown)>0){
        ## now check available paths
        next_path <- get_next_path(df_unknown[1,]$comb)
        if(!is.null(next_path)){
          next_comb <- get_xy_index(next_path, nrow(mat_phased))
        } else {
          next_comb <- NULL
        }      
      } else {
        next_comb <- NULL
      }      
    } else {
      next_comb <- NULL
    }
  } else {
    next_comb <- get_xy_index(shortest[1], nrow(mat_phased))
  }
  func_end()
  return(next_comb)
}
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom purrr set_names
#' @importFrom dplyr bind_rows left_join select filter
aggregate_phasing <- function(all_combs, df_gene, phasing_info){
  func_start()
  phasing_all_combs <- lapply(all_combs, function(mmp){
    comb_vec <- paste0("m",sort(get_xy_index(mmp, nrow(mat_phased))))
    comb <- comb_vec %>% paste(collapse="-")    
    df_gene_relcomb <- df_gene %>%
        filter(mut_id %in% comb_vec)
    min_tcn=min(df_gene_relcomb$tcn)
    nstatus <- mat_phased[mmp]
    status <- get_string_status(nstatus)      
    min_poss_wt_cp=calc_left_wt_copies(min_tcn,
                                         2,
                                         df_gene_relcomb$aff_cp[1],
                                         df_gene_relcomb$aff_cp[2])
    max_poss_wt_cp=calc_left_wt_copies(min_tcn,
                                         1,
                                         df_gene_relcomb$aff_cp[1],
                                         df_gene_relcomb$aff_cp[2])
    if(nstatus>0){
      ## status was defined
      ## calculate confidence
      used_combs <- mat_info[mmp] %>% 
        str_split("-") %>% 
        unlist() %>%
        as.numeric()
      extracted_combs <- phasing_info %>% filter(comb %in% used_combs) 
      conf <- extracted_combs %>%
        mutate(ep=ifelse(nstatus==2, 1-p_same, 1-p_diff)) %>%
        pull(ep) %>%
        prod()
      unplausible <- paste(as.numeric(extracted_combs$unplausible), collapse="-")
      subclonal <- paste(as.numeric(extracted_combs$subclonal), collapse="-")
      via <- as.character(mat_info[mmp])
      phasing=case_when(
        str_detect(mat_info[mmp], "h") ~ "haploblock",
        str_detect(mat_info[mmp], "s") ~ "imbalance",
        str_detect(mat_info[mmp], "-") ~ "indirect",
        TRUE ~ "direct"
      )
      wt_cp <- calc_left_wt_copies(min_tcn,
                                nstatus,
                                df_gene_relcomb$aff_cp[1],
                                df_gene_relcomb$aff_cp[2])
      score <- case_when(
        nstatus==2&wt_cp<0.5 ~ 2,
        TRUE ~ 1)
    } else {
      ## status could not be defined
      ## calculate maximum possible affected copies
      conf <- 0
      unplausible <- NA
      subclonal <- NA
      via <- NA
      wt_cp <- NA
      score <- case_when(
        min_poss_wt_cp>0.5 ~ 1,
        TRUE ~ 0
      )
      phasing <- case_when(
        score==1 ~ "exclusion",
        TRUE ~ NA
      )
    }
    phasing_status <- tibble(comb=comb, nstatus=nstatus, status=status, 
                             phasing=phasing, via=via, conf=conf, 
                             unplausible=unplausible,
                             subclonal=subclonal, wt_cp=wt_cp, 
                             min_poss_wt_cp=min_poss_wt_cp,
                             max_poss_wt_cp=max_poss_wt_cp,
                             score=score)
    return(phasing_status)
  }) %>%
    bind_rows()  %>%
    select(comb, nstatus, status, phasing, via, conf, unplausible, subclonal, 
           wt_cp, min_poss_wt_cp, max_poss_wt_cp, score)
  func_end()
  return(phasing_all_combs)
}
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom tibble tibble column_to_rownames as_tibble
#' @importFrom dplyr group_by tally mutate ungroup rowwise pull filter summarize left_join select
#' @importFrom stats chisq.test
classify_combination <- function(classified_reads, ref_class1, ref_class2, 
                                 purity, printLog){
  func_start()
  result <- . <- fac <- NULL
  all_possible_results <- c('both', 'mut1', 'mut2', 'none', 
                            "dev_var", 'read_in_read', 
                            'skipped')
  ## calculate expected counts according to aff copies for both cases:
  ## same and diff
  pseudo_count <- 0.0000001
  exp_diff <- tibble(fac=all_possible_results[1:3],
                     exp=c(0,1,1)+pseudo_count)
  exp_same <- tibble(fac=all_possible_results[1:3],
                     exp=c(1,0,0)+pseudo_count)
  number <- classified_reads %>%
    mutate(fac=factor(result, levels = all_possible_results)) %>%
    group_by(fac, .drop = FALSE) %>%
    tally() %>% 
    column_to_rownames(var='fac') %>%
    t() %>%
    as_tibble()
  none_raw <- number[['none']]
  append_loglist(print_tibble(number))
  ## calculate confidence from basecalls and mapping quality of read
  ## and aggregate them per classification result
  cr_conf <- classified_reads %>%
    rowwise() %>%
    mutate(pb1=ascii_to_dec(baseq1),
           pb2=ascii_to_dec(baseq2)) %>%
    ungroup() %>%
    mutate(
      mq1=10^(as.numeric(mapq1)/(-10)),
      mq2=10^(as.numeric(mapq2)/(-10)),
      ppos1=(1-pb1)*(1-mq1),
      ppos2=(1-pb2)*(1-mq2),
      p_result=ppos1*ppos2
    ) %>%
    mutate(fac=factor(result, levels = all_possible_results)) %>%
    group_by(fac, .drop = FALSE) %>%
    summarize(n=length(fac),
              prob_sum=sum(p_result))
  relevant_for_decision <- cr_conf %>%
    filter(fac %in% c("both","mut1", "mut2")) %>%
    select(fac, prob_sum)
  both <- relevant_for_decision %>%
    filter(fac=="both") %>%
    pull(prob_sum)
  mut1 <- relevant_for_decision %>%
    filter(fac=="mut1") %>%
    pull(prob_sum)
  mut2 <- relevant_for_decision %>%
    filter(fac=="mut2") %>%
    pull(prob_sum)
  append_loglist(print_tibble(relevant_for_decision))
  ## predefine null result if no evidence
  status <- "null"
  p <- 1
  xsq_same <- xsq_diff <- p_diff <- p_same <- v_same <- v_diff <- NA
  nstatus <- evidence <- certainty <- confidence <- conf_log <- 0
  subclonal <- unplausible <- FALSE
  ## get total number of relevant classifications
  sum_rel <- sum(relevant_for_decision$prob_sum)
  if(sum_rel==0){
    append_loglist("no evidence for any classification")
  } else if(both==0&mut1==0&ref_class2=="snp"){
    append_loglist("only mut2 detected, which is SNP")
  } else if(both==0&mut2==0&ref_class1=="snp"){
    append_loglist("only mut1 detected, which is SNP")
  } else {
    ## check for similarity of numbers in both, mut1 and mut2 by chi-squared
    suppressWarnings(
      sim_diff <- chisq.test(t(left_join(exp_diff, relevant_for_decision,
                                         by="fac") %>% 
                                 select(exp, prob_sum)))
    )
    suppressWarnings(
      sim_same <- chisq.test(t(left_join(exp_same, relevant_for_decision,
                                         by="fac") %>% 
                                 select(exp, prob_sum)))  
    )
    xsq_diff <- sim_diff$statistic[[1]]
    xsq_same <- sim_same$statistic[[1]]
    p_diff <- sim_diff$p.value
    p_same <- sim_same$p.value
    ## calculate cramers V... as we have always 2x3 table, m=1
    v_diff <- sqrt(xsq_diff/sum(sum_rel, exp_diff$exp))
    v_same <- sqrt(xsq_same/sum(sum_rel, exp_same$exp))
    if(p_diff>p_same){
      status <- "diff"
      nstatus <- 2
    } else {
      status <- "same"
      nstatus <- 1
    }  
  }
  ## check if evidence for subclonal sample is given  
  if(both>0&((mut1==0&mut2>0)|(mut1>0&mut2==0))){
    if(!"snp" %in% c(ref_class1, ref_class2)){
      warning("evidence fur subclonality found during phasing!")
      subclonal <- TRUE      
    }
  }
  ## check for unplausible case
  if(both>0&mut1>0&mut2>0){
    warning("unplausible phasing result!")
    unplausible <- TRUE
  }
  status_table <- tibble(
    both=both,
    mut1=mut1,  
    mut2=mut2,
    dev_var=number[['dev_var']],
    skipped=sum(as.numeric(number[['read_in_read']]), 
                   as.numeric(number[['skipped']])),
    status=status,
    nstatus=nstatus,
    xsq_same=xsq_same,
    xsq_diff=xsq_diff,
    p_same=p_same,
    p_diff=p_diff,
    v_same=v_same,
    v_diff=v_diff,
    none_raw=none_raw,
    DNA_rds=nrow(classified_reads %>% filter(origin=="DNA")),
    RNA_rds=nrow(classified_reads %>% filter(origin=="RNA")),
    subclonal=subclonal,
    unplausible=unplausible
  )
  func_end()
  return(status_table)
}  
#' @keywords internal
#' description follows
check_for_overlapping_reads <- function(bamDna, bamRna,
                                        ref_chr1, 
                                        ref_chr2, 
                                        ref_pos1, 
                                        ref_pos2){
  func_start()
  dna_bam <- prepare_raw_bam_file(bamDna, 
                                  ref_chr1, 
                                  ref_chr2, 
                                  ref_pos1, 
                                  ref_pos2) 
  if(length(dna_bam)!=0){
    dna_bam$origin <- "DNA"
  }
  if(!is.null(bamRna)){
    rna_bam <- prepare_raw_bam_file(bamRna, 
                                    ref_chr1, 
                                    ref_chr2, 
                                    ref_pos1, 
                                    ref_pos2)
    if(length(rna_bam)==0){
      rna_bam <- NULL
    } else {
      rna_bam$origin <- "RNA"
    }
  } else {
    rna_bam <- NULL
  }
  func_end()
  return(c(dna_bam, rna_bam))
}
#' @keywords internal
#' @importFrom magrittr %>%
#' @importFrom dplyr tibble bind_rows 
classify_reads <- function(ref_pos1,
                           ref_pos2,
                           ref_chr1,
                           ref_chr2,
                           ref_alt1,
                           ref_alt2,
                           ref_ref1,
                           ref_ref2,
                           ref_class1,
                           ref_class2, bamDna, bamRna){
  vm(as.character(sys.call()[1]),  1)
  bam <- check_for_overlapping_reads(bamDna,
                                     bamRna,
                                     ref_chr1,
                                     ref_chr2,
                                     ref_pos1,
                                     ref_pos2)
  if(!length(bam)==0){  
    classified_reads <- lapply(unique(bam$qname),
                               core_tool,
                               bam, 
                               ref_pos1,
                               ref_pos2,
                               ref_alt1,
                               ref_alt2,
                               ref_ref1,
                               ref_ref2,
                               ref_class1,
                               ref_class2) %>%
      bind_rows() 
  } else {
    classified_reads <- tibble()
  } 
  func_end()
  return(classified_reads)
}
#' @importFrom magrittr %>%
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
phase_combination <- function(mat_gene_relcomb, comb, bamDna, bamRna,  
                              geneDir, phasing_type, showReadDetail){
  func_start()
  append_loglist("Phasing:", comb)
  ## predefine empty output
  classified_main_comb <- tibble(comb=comb, status="null", nstatus=0, conf=0,
                                 phasing=comb)
  ## define relevant properties for mut1 and mut2 globally for the phasing process
  mut1 <- 2
  mut2 <- 1
  
  ref_pos1 <- as.numeric(mat_gene_relcomb[,"pos"][[mut1]])  
  ref_pos2 <- as.numeric(mat_gene_relcomb[,"pos"][[mut2]])
  ## must be character to work for annotations like: chrX
  ref_chr1 <- as.character(mat_gene_relcomb[,"chr"][[mut1]])  
  ref_chr2 <- as.character(mat_gene_relcomb[,"chr"][[mut2]])
  
  ref_alt1 <- as.character(mat_gene_relcomb[,"alt"][[mut1]])
  ref_alt2 <- as.character(mat_gene_relcomb[,"alt"][[mut2]]) 
  
  ref_ref1 <- as.character(mat_gene_relcomb[,"ref"][[mut1]])
  ref_ref2 <- as.character(mat_gene_relcomb[,"ref"][[mut2]])
  
  ref_class1 <- as.character(mat_gene_relcomb[,"class"][[mut1]])
  ref_class2 <- as.character(mat_gene_relcomb[,"class"][[mut2]])
  
  main_classified_reads <- classify_reads(ref_pos1,
                                          ref_pos2,
                                          ref_chr1,
                                          ref_chr2,
                                          ref_alt1,
                                          ref_alt2,
                                          ref_ref1,
                                          ref_ref2,
                                          ref_class1,
                                          ref_class2, 
                                          bamDna, 
                                          bamRna)
  append_loglist(nrow(main_classified_reads), 
                 "reads / read-pairs covering both positions")
  if(nrow(main_classified_reads)!=0){
    if(showReadDetail==TRUE){
      store_log(geneDir, main_classified_reads %>%
                  mutate(comb=comb),
                paste0("classified_reads_", comb,".tsv"))
    }
    classified_main_comb <- classify_combination(main_classified_reads,
                                                 ref_class1, ref_class2,
                                                 purity,
                                                 printLog
                                                 
    ) %>%  
      mutate(
        dist=mat_dist[comb],
        class_comb=paste(mat_gene_relcomb[,"class"][[2]], 
                         mat_gene_relcomb[,"class"][[1]], 
                         sep="-"),
        comb=comb,
        phasing=comb
      ) 
  }   
  func_end()
  return(classified_main_comb)
}
#' @importFrom magrittr %>%
append_matrices <- function(classified_main_comb, iterate=TRUE){
  func_start()
  if(nrow(classified_main_comb)>0){
    mat_old <- mat_phased
    lapply(seq(1,nrow(classified_main_comb)), function(C){
      mat_phased[classified_main_comb[C,]$comb] <<- 
        classified_main_comb[C,]$nstatus
      mat_info[classified_main_comb[C,]$comb] <<- 
        classified_main_comb[C,]$phasing
      return()
    })
    if(iterate==TRUE){
      something_changed <- sum(mat_old, na.rm = T)!=sum(mat_phased, na.rm = T)
      mat_new <- mat_phased
      mat_new[is.na(mat_new)] <- 0
      mat_old[is.na(mat_old)] <- 0
      changes <- which(!mat_new==mat_old) 
      i <- 1
      while(i<=length(changes)){
        m12 <- get_xy_index(changes[i], nrow(mat_new))
        status <- mat_new[changes[i]]
        ## get third mut/snp id that also has a non-null connection to one of the first ones 
        m1_same_conn <- c(m12[1], which(mat_new[m12[1],]==1),
                          which(mat_new[,m12[1]]==1)) %>% sort()
        m2_same_conn <- c(m12[2],
                          which(mat_new[m12[2],]==1),
                          which(mat_new[,m12[2]]==1)) %>% sort()
        if(sum(c(length(m1_same_conn), length(m2_same_conn)))>2){
          if(status==1){
            conns <- sort(unique(c(m1_same_conn, m2_same_conn)))
            for(x in 1:(length(conns)-1)){
              for(j in (x+1):length(conns)){
                p1 <- conns[x]
                p2 <- conns[j]
                notxf <- m12[which(!m12 %in% c(p1,p2))]
                if(length(notxf)>=1){
                  notx <- min(notxf)
                  if(is.na(mat_phased[p1,p2])|mat_phased[p1,p2]==0){
                    combined_info <- paste(
                      c(mat_info[
                        min(c(notx, p1)),
                        max(c(notx, p1))
                      ], mat_info[
                        min(c(notx, p2)),
                        max(c(notx, p2))                   
                      ]), 
                      collapse="-")
                    mat_info[p1,p2] <<- combined_info   
                  }
                } 
                if(is.na(mat_phased[p1,p2])|mat_phased[p1,p2]==0){
                  mat_phased[p1,p2] <<- 1
                }
              }
            }
          } else {
            for(x in m1_same_conn){
              for(j in m2_same_conn){
                sv <- sort(c(x,j))
                p1 <- sv[1]
                p2 <- sv[2]
                notxf <- m12[which(!m12 %in% c(x,j))]
                if(length(notxf)>=1){
                  notx <- min(notxf)
                  if(is.na(mat_phased[p1,p2])|mat_phased[p1,p2]==0){
                    mat_info[p1, p2] <<- paste(
                      c(mat_info[
                        min(c(notx, x)),
                        max(c(notx, x))
                      ], mat_info[
                        min(c(notx, j)),
                        max(c(notx, j))                   
                      ]), 
                      collapse="-")
                  }
                }
                if(is.na(mat_phased[p1, p2])|mat_phased[p1, p2]==0){
                  mat_phased[p1, p2] <<- 2
                }        
              }
            }
          }      
        }
        i <- i+1 
      } 
    }
  }
  func_end()
}
#' @importFrom tibble tibble
#' @importFrom dplyr  bind_rows
perform_direct_phasing <- function(mat_gene, bamDna, bamRna, purity, 
                                    printLog, showReadDetail, geneDir){
  func_start()
  phasing_type <- "direct"
  phasing_info <- tibble()
  unphased <- which(is.na(mat_phased))
  i <- 1
  while(i <= length(unphased)){
    ## as long as there are unphased combinations, try to phase
    comb <- unphased[i]
    relcombxy <- get_xy_index(comb, nrow(mat_phased))
    mat_gene_relcomb <- mat_gene[relcombxy,]
    
    classified_main_comb <- phase_combination(mat_gene_relcomb, comb, bamDna, bamRna, 
                                               geneDir, phasing_type, showReadDetail)
    phasing_info <- bind_rows(phasing_info,
                              classified_main_comb)
    i <- i+1
  }
  append_matrices(phasing_info)
  func_end()
  return(list(status=NULL,
              info=phasing_info,
              exit=tibble(phasing="direct", info="test")))
}
#' @keywords internal
#' description follows
#' @importFrom magrittr %>%
#' @importFrom stringr str_split str_detect
#' @importFrom dplyr case_when bind_rows mutate rowwise ungroup
#' @importFrom tibble as_tibble
phase <- function(df_gene,
                  somCna,  
                  bamDna, 
                  purity, 
                  sex, 
                  bamRna=NULL, 
                  haploBlocks=NULL, 
                  vcf=NULL,
                  distCutOff=5000, 
                  printLog=FALSE, 
                  logDir=NULL, 
                  showReadDetail=FALSE, 
                  snpQualityCutOff=1, 
                  phasingMode="full", 
                  verbose=FALSE){
  vm(as.character(sys.call()[1]),  1)
  #haploblock_phasing <- NULL
  GENE <- unique(df_gene$gene)
  if(!is.null(logDir)){
    geneDir <- file.path(logDir, GENE)
    dir.create(geneDir)
  } else {
    geneDir <- NULL
  }
  ## (1): define all combinations of variants to be phased
  create_phasing_matrices(df_gene$mut_id, df_gene$pos, distCutOff)
  solved_master_combs <<- list()
  mat_gene <- as.matrix(df_gene)
  rownames(mat_gene) <- mat_gene[,"mut_id"]
  unphased <- unphased_main()
  append_loglist(length(unphased), "main combinations to phase,", 
                 abs(length(unphased)-length(mat_phased[upper.tri(mat_phased)])), 
                 "are over distCutOff")
  ## (2): perform direct phasing between variants (read-based)
  direct_phasing <- perform_direct_phasing(mat_gene, bamDna, bamRna, 
                                           purity,  printLog, 
                                           showReadDetail, geneDir)
  phasing_info <- direct_phasing$info
  if(length(unknown_main())>0&!is.null(vcf)){
    append_loglist("unphased combinations left --> Initialize SNP phasing")
    ## missing combinations in main muts --> start secondary phasing approaches
    lsnps <- load_snps(df_gene, vcf, haploBlocks,  distCutOff,  somCna, snpQualityCutOff)
    store_log(geneDir, lsnps, "lsnps.tsv")
    if(!is.null(lsnps)){
      snps <- lsnps %>%
        ## filter low quality
        rowwise() %>%
        mutate(
          aff_cp=ifelse(is.na(af)|is.na(tcn),
                        NA,
                        aff_germ_copies(seqnames, af, tcn, purity, sex))
        ) %>%
        ungroup() %>%
        mutate(
          gt_cna_min=str_split(gt_cna, ":") %>% unlist() %>% min(),
          gt_cna_max=str_split(gt_cna, ":") %>% unlist() %>% max(),
          gt_seg=case_when(
            all_imb==FALSE|is.na(aff_cp) ~ NA,
            round(aff_cp) == gt_cna_min ~ "0|1",
            round(aff_cp) == gt_cna_max ~ "1|0",
          ),
          gt_final=case_when(
            gt=="1|1"|gt=="0|0" ~ NA,
            !is.na(hap_id)&str_detect(gt, "\\|") ~ gt,
            !is.na(gt_seg) ~ gt_seg,
            TRUE ~ NA
          ),
          block_final=case_when(
            gt=="1|1"|gt=="0|0" ~ NA,
            !is.na(hap_id)&str_detect(gt, "\\|") ~ paste0("h", hap_id),
            !is.na(gt_seg) ~ paste0("s", seg_id),
            TRUE ~ NA

          )
        )
      append_phasing_matrices(snps$mut_id, snps$start, distCutOff)
      add_snps_to_matrices(snps)
      ## build main daatframe of muts and snps and annotate haploblocks and segments of allelic imbalance
      df <- bind_rows(
        df_gene %>% select(chr, pos, ref, alt, class, af, mut_id),
        snps %>% select(chr=seqnames, pos=start, ref=REF, alt=ALT, af, mut_id, gt=gt_final) %>% mutate(class="snp")
      )
      to_phase <- prioritize_combination()
      i <- 1
      while(!is.null(to_phase)&length(unknown_main())>0){
        append_loglist("Phasing combination:", paste(paste0(to_phase), collapse="-"))

        mat_gene_relcomb <- df[to_phase,] %>% as.matrix()
        comb <- get_single_index(to_phase[1],
                                 to_phase[2],
                                 nrow(mat_phased))
        classified_main_comb <- phase_combination(mat_gene_relcomb, comb,
                                                  bamDna, bamRna,  
                                                  geneDir, comb, showReadDetail)
        append_matrices(classified_main_comb)
        phasing_info <- bind_rows(phasing_info,
                                  classified_main_comb)
        to_phase <- prioritize_combination()
        i <- i +1
      }
      store_log(geneDir, df, "df_snps_muts.tsv")
      store_log(geneDir, phasing_info, "all_indirect_phasing_combinations.tsv")
    }## no snps found
  }
  ## secondary phasing done
  append_loglist("finalizing phasing results")
  ## reconstruct phasing results
  main_pos <- get_main_mut_pos()
  uppertri <- which(upper.tri(mat_phased))
  all_combs <- intersect(main_pos, uppertri)
  phasing_all_combs <- aggregate_phasing(all_combs, df_gene, phasing_info)
 if(nrow(phasing_info)==0){
   phasing_info_export <- tibble()
 } else {
    phasing_info_export <- phasing_info %>%
      rowwise() %>%
      mutate(ncomb=comb,
             gene=GENE,
             comb=paste(
               sort(
                 factor(
                  c(colnames(mat_phased)[get_xy_index(ncomb, 
                                                      nrow(mat_phased))["x"]], 
                     rownames(mat_phased)[get_xy_index(ncomb, 
                                                       nrow(mat_phased))["y"]]),
                  levels=c(paste0("m", seq(1,100,1)),
                           paste0("s", seq(1,2000,1)))
                 )
               ),
               collapse="-"
               ))   
  }
  store_log(geneDir, as_tibble(mat_phased), "mat_phased.tsv")
  store_log(geneDir, as_tibble(mat_info), "mat_info.tsv")
  store_log(geneDir, phasing_info, "all_phasing_combinations.tsv")
  func_end()
  return(list(phasing_all_combs, phasing_info_export, mat_phased, mat_info))
  
}
#' @keywords internal
calc_left_wt_copies <- function(mtcn, nstatus, aff_copies1, aff_copies2){
  left_wt_copies <- ifelse(nstatus==2,
                           mtcn-sum(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)),
                           mtcn-max(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)))
  return(as.numeric(left_wt_copies))
}
#' @keywords internal
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_detect str_replace_all
#' @importFrom dplyr mutate
#' @importFrom tibble as_tibble
make_both_annotations <- function(region_to_load_in){
  sec <- region_to_load_in %>%
    as_tibble()
  if(str_detect(sec$seqnames[1],"chr")){
    prim <- sec %>%
      mutate(seqnames=str_replace_all(seqnames, "chr", "")) %>%
      GRanges()
  } else {
    prim <- sec %>%
      mutate(seqnames=paste0("chr",seqnames)) %>%
      GRanges()
  }
  fin <- c(region_to_load_in, prim)
  return(fin)
}
#' @importFrom GenomicRanges GRanges
#' @importFrom DelayedArray rowRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom stringr str_detect str_replace str_replace_all
#' @importFrom dplyr mutate pull nth 
#' @importFrom tibble as_tibble
#' @importFrom purrr compact
#' @importFrom VariantAnnotation readVcf geno info
#' @importFrom Rsamtools seqnamesTabix TabixFile
loadVcf <- function(vcf_in, chrom, region_to_load_in,  somCna,which="all",
                    colname_gt="GT", colname_af="AF", colname_dp4="DP4", 
                    dkfz=FALSE){
  func_start()
  ## first check which format input vcf has
  region_to_load <- region_to_load_in
  chr_anno <- as_tibble(region_to_load_in) %>%
    pull(seqnames) %>% .[1] %>% str_detect("chr")
  sec_chrom <- ifelse(chr_anno, str_replace(chrom, "chr", ""), 
                      paste0("chr", chrom))
  gr_list <- lapply(vcf_in, function(VCF){
    if(is(VCF, "TabixFile")){
      if(chrom %in% Rsamtools::seqnamesTabix(VCF)){
        loadedVcf <- 
          readVcf(VCF, 
                  param=region_to_load)
        combVcf <- rowRanges(loadedVcf) 
        if(length(combVcf)>0){
          combVcf$ALT <- unlist(lapply(combVcf$ALT, 
                                       function(x){as.character(x[[1]][1])}))
          gt <- VariantAnnotation::geno(loadedVcf)[[colname_gt]][,1] %>% 
            as.character()
          af <- VariantAnnotation::info(loadedVcf)[[colname_af]]
          dp4 <- VariantAnnotation::info(loadedVcf)[[colname_dp4]] 
          combVcf$gt <- gt
          combVcf$dp4 <- dp4
          combVcf$af <- af
          vcf_info <- "variants_detected"
        } else {
          combVcf <- NULL
          vcf_info <- "no snps in roi"
        }
      } else {
        combVcf <- NULL
        vcf_info <- "vcf does not cover chromosome"
      }
    } else {
      if(str_detect(VCF, paste0("chr", 
                                unlist(str_replace(chrom, "chr", "")), 
                                "[^0-9]"))){
        message("loading")
        loadedVcf <- VariantAnnotation::readVcf(VCF)
        gt <- VariantAnnotation::geno(loadedVcf)[[colname_gt]] %>% 
          as.character()
        af <- VariantAnnotation::info(loadedVcf)[[colname_af]] 
        dp4 <- VariantAnnotation::info(loadedVcf)[[colname_dp4]]
        rangesVcf <- rowRanges(loadedVcf) 
        rangesVcf$gt <- gt
        rangesVcf$af <- af
        rangesVcf$dp4 <- dp4
        combVcf <-  subsetByOverlaps(rangesVcf, region_to_load) 
        if(length(combVcf)>0){
          combVcf$ALT <- unlist(lapply(combVcf$ALT, 
                                       function(x){as.character(x[[1]][1])}))
          combVcf <- as_tibble(combVcf) %>%
            mutate(seqnames=case_when(
              chr_anno==TRUE&!str_detect(seqnames, "chr") ~ paste0("chr",seqnames),
              chr_anno==FALSE&str_detect(seqnames, "chr") ~ str_replace_all(seqnames, "chr", ""),
              TRUE ~ seqnames
              
            )) %>%
            GRanges()
          vcf_info <- "variants_detected"
        } else {
          combVcf <- NULL
          vcf_info <- "no snps in roi"
        }
      } else {
        combVcf <- NULL
        vcf_info <- "no vcf file covering roi"
      }
    }
    return(list(combVcf, vcf_info))
  }) 
  final_vcf <- lapply(gr_list, nth, 1) %>%
    compact() %>%
    Reduce(function(x,y)c(x,y),.)
  if(which=="HZ"){
    filtered_vcf <- final_vcf[which(final_vcf$gt %in% c("1/0", "0/1", "1|0", "0|1"))]
  } else {
    filtered_vcf <- final_vcf
  }
  func_end()
  return(filtered_vcf)
}
#' @importFrom dplyr case_when
get_string_status <- function(nstatus){
  case_when(nstatus==2 ~ "diff",
            nstatus==1 ~ "same",
            TRUE ~ "null"
  )
}
#' @importFrom IRanges mergeByOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr filter select mutate left_join
#' @importFrom tibble as_tibble
load_snps <- function(df_gene, vcf, haploBlocks, distCutOff,  somCna, 
                      snpQualityCutOff, which="all"){
  func_start()
  region_to_load <- paste0(unique(df_gene$chr), ":", min(df_gene$pos)-distCutOff,
                           "-", max(df_gene$pos)+distCutOff) %>%
    GRanges()  %>%
    split_genomic_range(.,df_gene$pos)
  loaded_vcf_hz <- loadVcf(vcf, unique(df_gene$chr), region_to_load,  
                            somCna)
  if(length(loaded_vcf_hz)>0){
    fsnps <- loaded_vcf_hz %>% as_tibble() %>%
      filter(!is.na(QUAL)&QUAL>snpQualityCutOff) %>%
      filter(gt %in% c("1|0", "0|1", "1/0", "0/1"))
    if(nrow(fsnps)>0){
      snps <- fsnps %>%
        mutate(mut_id=paste0("s", c(1:nrow(.))+nrow(df_gene)))
      append_loglist(nrow(snps), "SNPs detected in distCutOff range")
      if(!"af" %in% names(snps)){
        if("dp4" %in% names(snps)){
          snps$af <- lapply(seq(1,nrow(snps)), function(i){
            sum(snps$dp4[[i]][c(3,4)])/sum(snps$dp4[[i]])
          }) %>% as.numeric()  
        } 
      }
      gr_snps <- GRanges(snps %>% select(1:3, mut_id))
      if(!is.null(haploBlocks)){
        merged_haploblocks <- mergeByOverlaps(gr_snps, haploBlocks) %>%
          as_tibble() %>%
          select(all_of(names(.) %>% .[which(!str_detect(.,"\\."))]))    
      } else {
        merged_haploblocks <- snps %>% select(mut_id) %>%
          mutate(hap_id=NA)
      }
      merged_somCna <- mergeByOverlaps(gr_snps, somCna) %>%
        as_tibble() %>%
        select(all_of(names(.) %>% .[which(!str_detect(.,"\\."))]))
      annotated_snps <- left_join(
        snps,
        merged_haploblocks,
        by="mut_id"
      ) %>%
        left_join(
          .,
          merged_somCna,
          by="mut_id"
        )
    } else {
      ## no high quality and heterozygous snps
      annotated_snps <- NULL
    }
  } else {
    ## no snps detected
    annotated_snps <- NULL
  }
  func_end()
  return(annotated_snps)
}

get_genotype <- function(gt, status){
  if(status=="same"){
    return(gt)
  } else if(gt=="1|0"){
    return("0|1")
  } else if(gt=="0|1"){
    return("1|0")
  } 
}
aggregate_probs <- function(ps){
  if(length(ps)>1){
    np <- ps[1]
    for (n in seq(1, length(ps)-1, 1)){
      np <- np+(1-np)*ps[n+1]
    }
    return(np)    
  } else {
    return(ps) 
  }
}
#' @importFrom stringr str_sub
#' @importFrom stats na.omit
#' @importFrom dplyr pull all_of
extract_subseq <- function(ref_pos, element_mut, parsed_read, 
                           length_indel, string){
  pos_in_seq <- ref_pos-element_mut$map_start+1
  subseq_in_element <- str_sub(element_mut[[string]], start=pos_in_seq, 
                               end=-1)
  subseq_in_following_elements <- 
    parsed_read[
      which(as.numeric(parsed_read$id)>as.numeric(element_mut$id)),] %>%
    pull(all_of(string)) %>%
    na.omit() %>%
    paste(collapse="")
  full_remaining <- paste0(subseq_in_element, subseq_in_following_elements)
  return(str_sub(full_remaining, start=1, end=length_indel))
}
#' @importFrom stringr str_sub
#' @importFrom tibble tibble
extract_snv <- function(ref_pos, element_mut){
  ## get cigar element in which reference position in located
  pos_in_seq <- ref_pos-element_mut$map_start+1
  base <- 
    str_sub(element_mut$seq,
            start=pos_in_seq,
            end=pos_in_seq) 
  qual <- 
    str_sub(element_mut$qual,
            start=pos_in_seq,
            end=pos_in_seq)  
  return(tibble(base=base, qual=qual, info="mapped", len_indel=NA,
                exp_indel=NA))
}
#' @importFrom stringr str_sub
#' @importFrom tibble tibble
extract_insertion <- function(ref_pos, parsed_read, element_mut, length_ins){
  ## insertion are always indicated by the I type. The reference position tells
  ## the position in front of the I segment, which means we are looking for
  ## the subsequent one type
  expected_indel_detected <- FALSE
  len <- NA
  if(as.numeric(element_mut$id)==nrow(parsed_read)){
    if(ref_pos==element_mut$map_end){
      info <- 
        "no insertion detected: ref pos is last mapped base of read" 
    } else {
      info <- 
        "no insertion detected: ref pos mapped in last element of parsed read"      
    }
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    #qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  } else if("I" %in% parsed_read$type){
    ## insertion detected... check if it is the matching one
    if(ref_pos==element_mut$map_end){
      ## if insertion is at refrence position, the subsequent position is an a
      ## new elemnt labelled as I
      next_element <- parsed_read[as.numeric(element_mut$id)+1,]
      if(next_element$type=="I"){
        ## insertion detected at position.. extract inserted sequence
        base <- paste0(
          str_sub(element_mut$seq, start=-1, end=-1),
          next_element$seq
        )
        #  qual <- paste0(
        #   str_sub(element_mut$qual, start=-1, end=-1),
        #  next_element$qual
        #)
        info <- "insertion detected"
        len <- next_element$width
        expected_indel_detected <- TRUE
      } else {
        info <- paste(next_element$type, "detected - wrong class")
        base <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
                               "seq")
        #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_ins, 
        #                      "qual")
      }
    } else {
      info <- "no insertion detected: ref pos not mapped to end of element"
      base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
      # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
    }
  } else {
    ## no inserttion detected.. extract base at position to get base quality
    info <- "no insertion detected: no I in cigar"
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  }
  ## qual is set to null, as we only take the mapping quality for indels
  return(tibble(base=base, qual=0, info=info, len_indel=len, 
                exp_indel=expected_indel_detected))
}
#' @importFrom tibble tibble
extract_deletion <- function(ref_pos, parsed_read, element_mut, length_del){
  expected_indel_detected <- FALSE
  len <- NA
  if(as.numeric(element_mut$id)==nrow(parsed_read)){
    if(ref_pos==element_mut$map_end){
      info <- 
        "no deletion detected: ref pos is last mapped base of read" 
    } else {
      info <- 
        "no delection detected: ref pos mapped in last element of parsed read"      
    }
    base <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "seq")
    # qual <- extract_subseq(ref_pos, element_mut, parsed_read, 1, "qual")
  } else if("D" %in% parsed_read$type){
    ## insertion detected... check if it is the matching one
    if(ref_pos==element_mut$map_end){
      ## if insertion is at refrence position, the subsequent position is an a
      ## new elemnt labelled as D
      next_element <- parsed_read[as.numeric(element_mut$id)+1,]
      if(next_element$type=="D"){
        expected_indel_detected <- TRUE
        info <- "deletion detected"
        len <- next_element$width
      } else {
        info <- "no deletion detected: subsequent element not I"
      }
    } else {
      info <- "no insertion detected: ref pos not mapped to end of element"
    }      
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "seq")
    #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
    #                      "qual")
  } else {
    ## no inserttion detected.. extract base at position to get base quality
    info <- "no insertion detected: no I in cigar"
    base <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
                           "seq")
    #qual <- extract_subseq(ref_pos, element_mut, parsed_read, length_del, 
    #                      "qual")
  }
  ## qual is set to null, as we only take the mapping quality for indels
  return(tibble(base=base, qual=0, info=info, len_indel=len,
                exp_indel=expected_indel_detected))
}
#' @importFrom dplyr between
cigar_element <- function(parsed_read, ref_pos){
  #func_start()
  cigel <- parsed_read[which(between(rep(ref_pos, nrow(parsed_read)), 
                                     parsed_read$map_start, 
                                     parsed_read$map_end)),] %>%
    .[1,]
  #func_end()
  return(cigel)
}
#' @importFrom tibble tibble
extract_base_at_refpos <- function(parsed_read, ref_pos, class, ref_alt, 
                                   ref_ref){
  #func_start()
  ## extract elemnt containing reference position from parsed read
  element_mut <- cigar_element(parsed_read, ref_pos)
  ## if position is inside N or D type, the read does not cover the position
  ## N means for example skipped exon in RNA
  ## D wouldmean the position is deleted
  ## real deletions are still mapped in the M part
  if(element_mut$type %in% c("N", "D")){
    base_info <- tibble(base=NA, qual=NA, info="position skipped/deleted", 
                        len_indel=NA, exp_indel=NA)
  } else if(class %in% c("snv", "snp")){
    base_info <- extract_snv(ref_pos, element_mut)
  } else if(class=="ins"){
    base_info <- extract_insertion(ref_pos, parsed_read, element_mut, 
                                   nchar(ref_alt)-nchar(ref_ref))
  } else if(class=="del"){
    base_info <- extract_deletion(ref_pos, parsed_read, element_mut, 
                                  1+abs(nchar(ref_ref)-nchar(ref_alt)))
  } else {
    stop("class of variant needs to be provided")
  }
  #func_end()
  return(base_info %>% mutate(class=class, mapq=element_mut$mapq))
}
#' @importFrom stringr str_sub
#' @importFrom dplyr bind_rows mutate na_if
#' @importFrom GenomicAlignments start end width
parse_cigar <- function(bam, qname){
  #func_start()
  paired_reads <- bam[which(bam$qname==qname)] 
  read_start <- GenomicAlignments::start(paired_reads)
  seq <- as.character(paired_reads$seq)
  cigar <- paired_reads$cigar
  qual <- as.character(paired_reads$qual)
  mate <- c(1,2)
  ## parse cigar string according to query
  cigq <- GenomicAlignments::cigarRangesAlongQuerySpace(cigar,
                                                        with.ops = T) 
  ## parse cigar string according to reference
  cigr <- GenomicAlignments::cigarRangesAlongReferenceSpace(cigar,
                                                            with.ops = F) 
  raw_cigs <- lapply(mate, function(i){
    cigq_start <- GenomicAlignments::start(cigq[[i]])
    cigq_end <- GenomicAlignments::end(cigq[[i]])
    cigq_len <- length(cigq_start)
    cigr_start <- GenomicAlignments::start(cigr[[i]])
    cigr_end <- GenomicAlignments::end(cigr[[i]])
    seqq <- str_sub(seq[i],
                    start=cigq_start,
                    end=cigq_end) %>% na_if("")
    qualq <- str_sub(qual[i],
                     start=cigq_start,
                     end=cigq_end) %>% na_if("")
    width=pmax(GenomicAlignments::width(cigq[[i]]), 
               GenomicAlignments::width(cigr[[i]]))
    map_start=read_start[i]+cigr_end-width
    map_end=read_start[i]+cigr_end-1
    raw_cigs_new <- data.frame(
      width=width,
      type=names(cigq[[i]]),
      seq=seqq,
      qual=qualq,
      map_start=map_start,
      map_end=map_end,
      start=cigr_start,
      end=cigr_end,
      mate=mate[i],
      mapq=paired_reads$mapq[i],
      origin=paired_reads$origin[i]
    )
  }) %>%
    bind_rows() %>%
    mutate(id=c(1:nrow(.)))
  #func_end()
  return(raw_cigs)
}
#' @importFrom dplyr case_when
evaluate_base <- function(base_info, ref_alt, ref_ref){
  if(base_info$class %in% c("snv", "snp")){
    detected <- case_when(
      base_info$base==ref_alt ~ 1,
      base_info$base==ref_ref ~ 0,
      TRUE ~ -2
    )
  } else if(base_info$class=="ins"){
    detected <- case_when(
      base_info$base==ref_alt ~ 1,
      base_info$base==ref_ref ~ 0,
      TRUE ~ -2
    )
  } else {
    detected <- case_when(
      base_info$len_indel==abs(nchar(ref_alt)-nchar(ref_ref)) ~ 1,
      is.na(base_info$len_indel) ~ 0,
      TRUE ~ -2
    )
  }
  return(detected)
}
#' @keywords internal
#' description follows
#' @importFrom dplyr case_when 
core_tool <- function(qname, bam,
                      ref_pos1, ref_pos2,
                      ref_alt1, ref_alt2,
                      ref_ref1, ref_ref2,
                      ref_class1, ref_class2, version="old"){
  #vm("core_tool",  1)
  . <- NULL
  ## parse read according to cigar string
  parsed_read <- parse_cigar(bam, qname)
  ## extract base at reference position
  base_info1 <- extract_base_at_refpos(parsed_read, ref_pos1, ref_class1, 
                                       ref_alt1, ref_ref1)
  base_info2 <- extract_base_at_refpos(parsed_read, ref_pos2, ref_class2, 
                                       ref_alt2, ref_ref2)
  ## assign final status to read
  if(is.na(base_info1$base)|is.na(base_info2$base)){
    final_assignment <- "skipped"
  } else {
    mut1_in_read <- evaluate_base(base_info1, ref_alt1, ref_ref1)
    mut2_in_read <- evaluate_base(base_info2, ref_alt2, ref_ref2)
    final_assignment <- case_when(
      sum(mut1_in_read, mut2_in_read)==2 ~ "both",
      sum(mut1_in_read, mut2_in_read)==0 ~ "none",
      mut1_in_read==1 ~ "mut1",
      mut2_in_read==1 ~ "mut2",
      TRUE ~ "dev_var"
    )    
  }
  return(c(qname=qname, result=final_assignment,
           origin=unique(parsed_read$origin),
           baseq1=base_info1$qual,
           mapq1=base_info1$mapq,
           baseq2=base_info2$qual,
           mapq2=base_info2$mapq
  )
  )
}