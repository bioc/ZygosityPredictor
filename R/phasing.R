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
      #print(snps_in_hap)
      for (i in 1:(nrow(snps_in_hap)-1)) {
        result_matrix <- matrix(0, nrow = nrow(snps_in_hap)-i, ncol = 2)
        for (j in (i+1):nrow(snps_in_hap)) {
          #message(i, " - ", j)
          comparison_result <- ifelse(snps_in_hap$gt[i] == snps_in_hap$gt[j], 1, 2)
          # Store the results in the result_matrix
          result_matrix[nrow(snps_in_hap)-j+1, 1] <- get_single_index(snps_in_hap$id[j], snps_in_hap$id[i], nrow(mat_phased))
          result_matrix[nrow(snps_in_hap)-j+1, 2] <- comparison_result
          
        }
        res_list[[i]] <- result_matrix
      }
      
      to_add <- Reduce(function(x,y)rbind(x,y),res_list) #%>% t() %>%
      colnames(to_add) <- c("comb", "nstatus") 
      new <- to_add %>%
        as_tibble() %>%
        mutate(conf=5, phasing=HAP_ID, hap_id=HAP_ID)     
    } else {
      ## if haploblock has only one snp, no concludions can be drawn
      new <- tibble()
    }
    #%>
    
    return(new)
  })  %>% bind_rows()
  #print(22222222222222222222222)
  #print(mat_phased)
  #print(known_combs)
  append_matrices(known_combs, iterate=FALSE)
  #print(mat_phased)
  func_end()
}
create_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  mat_dist <<- make_dist_matrix(all_pos, all_variants, distCutOff) 
  #if(sum(mat_dist[seq(1,length(df_gene$mut_id)),])>0){
  main_muts <<- all_variants
  main_pos <<- all_pos
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA    
  #mat_conf <<- mat_phased
  mat_info <<- mat_phased
  #} else {
  ## no snps closer than distCutoff to main muts
  #}
}
append_phasing_matrices <- function(all_variants, all_pos, distCutOff){
  func_start()
  mat_dist <<- make_dist_matrix(c(main_pos, all_pos), 
                                c(main_muts, all_variants), 
                                distCutOff) 
  mat_phased_main <- mat_phased
  #mat_conf_main <- mat_conf
  mat_info_main <- mat_info
  
  
  
  mat_phased <<- mat_dist
  mat_phased[mat_phased!=0] <<- NA 
  #mat_conf <<- mat_phased
  mat_info <<- mat_phased
  #print(mat_phased)
  
  nm <- length(main_muts)
  #print(mat_phased[nm, nm])
  #print(mat_phased_main)
  mat_phased[c(1:nm), c(1:nm)] <<- mat_phased_main
  #mat_conf[c(1:nm), c(1:nm)] <<- mat_conf_main
  mat_info[c(1:nm), c(1:nm)] <<- mat_info_main
  func_end()
}
get_main_mut_conns <- function(){
  lapply(main_muts, function(M){
    length(c(which(mat_phased[M,]!=0|is.na(mat_phased[M,])),
             which(mat_phased[,M]!=0|is.na(mat_phased[,M]))
    ))
  }) %>% unlist() %>%
    set_names(main_muts)
}
get_next_path <- function(comb, distCutOff){
  #if(sum(is.na(mat_phased))>0){
  func_start()
      
      mains <- as.character(sort(get_xy_index(as.numeric(comb), nrow(mat_phased))))
  
      # main_mut_conns <- get_main_mut_conns()
      # nn_mmc <- main_mut_conns[which(main_mut_conns!=0)]
   # if(length(nn_mmc)>1){
      ldc <- 1000
      #while(length(all_paths)==0&ldc<distCutOff){
        
        #print(ldc)
        open_conns <- which((is.na(mat_phased)|mat_phased>0)&mat_dist<ldc) %>%
          get_xy_index(.,nrow(mat_phased))
       #print(open_conns)
        if(length(dim(open_conns))!=2){
          connections <- as_tibble(t(open_conns)) %>%
            set_names(c("mut_id1", "mut_id2"))
        } else {
          colnames(open_conns) <- c("mut_id1", "mut_id2")
          connections <- as_tibble(open_conns)
        }
  if(length(intersect(as.numeric(mains), c(connections$mut_id1, connections$mut_id2)))==2){
     #print(connections)
  
       #print(connections$mut_id1)
       #print(connections$mut_id2)
       #print(mains)
        graph <- graph_from_data_frame(connections, directed = FALSE)
       #print(graph)
        #shortest_path <-
        #nnn_mmc <- names(nn_mmc) %>% str_match("\\d") %>% unlist()
        #print(nnn_mmc)
        all_paths <- igraph::all_shortest_paths(graph, from = mains[1], to = mains[2], mode = "all")$res  
       #ldc <- ldc*2      
      #}


      if(length(all_paths)>0){
        poss <- list()
        for (PATH in all_paths){
          P <- as.numeric(names(PATH))
          res <- c()
          for(i in c(1:(length(P)-1))){
            #print(i)
            # for(j in c((i+1):length(P))){
            #message(i, "-", i+1)
            ##print(j)
            #print(P[i])
            #print(P[j])
            mv <- sort(c(P[i], P[i+1]))
            #print(mv)
            conn <- get_single_index(mv[2], mv[1], nrow(mat_phased))
            #print(conn)
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
prioritize_combination <- function(){
  #print(1)
  func_start()
  unrel_rows <- rownames(mat_phased) %>% .[which(!. %in% colnames(mat_phased))]
  mat_tmp <- mat_dist
  mat_tmp[unrel_rows,] <- NA
  mat_tmp[,unrel_rows] <- NA
  #print(main_muts)
  #print(mat_phased)
  mat_tmp[c((length(main_muts)+1):nrow(mat_phased)),] <- NA
  #print(2)
  unphased <- which(is.na(mat_phased))
  relevant_dists <- mat_tmp[unphased] 
  relevant_dists[which(relevant_dists==0)] <- NA
  #print(3)
  shortest <- unphased[which(relevant_dists==min(relevant_dists, na.rm=T))] 
  ## pick first one, if two have the same distance
  #print(shortest)
  #print(mat_phased)
  if(length(shortest)==0){
    #stop("stop_here")
    main_mut_conns <- get_main_mut_conns()
    df_mmc <- tibble(m=names(main_mut_conns),
                     n=main_mut_conns)
    uc <- unknown_main()
    if(length(uc)>0){
      #print(uc)
      #print(mat_phased)
      #print(df_mmc)
      df_unknown_pre <- lapply(uc, function(C){
        sort(get_xy_index(C, nrow(mat_phased))) %>%
          paste0("m",.) %>%
          set_names(c("m1", "m2")) %>%
          c(.,comb=C)
      }) %>%
        bind_rows() #%>%
      #print(df_unknown_pre)
      df_unknown <- df_unknown_pre %>%
        #mutate_at(.vars = c("m1", "m2"), .funs = paste0) %>%
        left_join(df_mmc %>% select(m, n1=n), by=c("m1"="m")) %>%
        left_join(df_mmc %>% select(m, n2=n), by=c("m2"="m")) %>%
        filter(n1>0&n2>0)
     
     # print(df_unknown)
      
      if(nrow(df_unknown)>0){
        ## now check available paths
        next_path <- get_next_path(df_unknown[1,]$comb)
        
        if(!is.null(next_path)){
          #print(next_path)
          next_comb <- get_xy_index(next_path, nrow(mat_phased))
        } else {
          next_comb <- NULL
        }      
      } else {
        next_comb <- NULL
      }      
    }else {
      next_comb <- NULL
    }

  } else {
    next_comb <- get_xy_index(shortest[1], nrow(mat_phased))
  }
  func_end()
  return(next_comb)
}
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
      
      #print(used_combs)
      conf <- extracted_combs %>%
        mutate(ep=ifelse(nstatus==2, 1-p_same, 1-p_diff)) %>%
        pull(ep) %>%
        prod()
      
      unplausible <- paste(as.numeric(extracted_combs$unplausible), collapse="-")
      subclonal <- paste(as.numeric(extracted_combs$subclonal), collapse="-")
      #print(11111111111111111111)
      #print(conf)
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

      
     # min_poss_wt_cp <- NA
    #  max_poss_wt_cp <- NA
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
                             phasing=phasing, via=via, conf=conf, unplausible=unplausible,
                             subclonal=subclonal, wt_cp=wt_cp, min_poss_wt_cp=min_poss_wt_cp,
                             max_poss_wt_cp=max_poss_wt_cp,
                             score=score)
    return(phasing_status)
  }) %>%
    bind_rows()  %>%
    select(comb, nstatus, status, phasing, via, conf, unplausible, subclonal, wt_cp, min_poss_wt_cp, max_poss_wt_cp, score)
  func_end()
  return(phasing_all_combs)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr as_tibble group_by mutate tibble tally 
#' @importFrom purrr set_names
classify_combination <- function(classified_reads, ref_class1, ref_class2, 
                                 purity, printLog, 
                                 verbose=FALSE){
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
    sim_diff <- chisq.test(t(left_join(exp_diff, relevant_for_decision,
                                       by="fac") %>% 
                               select(exp, prob_sum)))
    sim_same <- chisq.test(t(left_join(exp_same, relevant_for_decision,
                                       by="fac") %>% 
                               select(exp, prob_sum)))     
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
  func_end(verbose)
  return(status_table)
}  
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr tibble 
classify_reads <- function(ref_pos1,
                           ref_pos2,
                           ref_chr1,
                           ref_chr2,
                           ref_alt1,
                           ref_alt2,
                           ref_ref1,
                           ref_ref2,
                           ref_class1,
                           ref_class2, bamDna, bamRna, verbose){
  vm(as.character(sys.call()[1]), verbose, 1)
  #dist <- abs(as.numeric(line[['pos1']])-as.numeric(line[['pos2']]))
  
  # ref_pos1 <- as.numeric(mat_gene_relcomb[,"pos"][[2]])  
  # ref_pos2 <- as.numeric(mat_gene_relcomb[,"pos"][[1]])
  # 
  # ref_chr1 <- as.numeric(mat_gene_relcomb[,"chr"][[2]])  
  # ref_chr2 <- as.numeric(mat_gene_relcomb[,"chr"][[1]])
  # 
  # ref_alt1 <- as.character(mat_gene_relcomb[,"alt"][[2]])
  # ref_alt2 <- as.character(mat_gene_relcomb[,"alt"][[1]]) 
  # 
  # ref_ref1 <- as.character(mat_gene_relcomb[,"ref"][[2]])
  # ref_ref2 <- as.character(mat_gene_relcomb[,"ref"][[1]])
  # 
  # ref_class1 <- as.character(mat_gene_relcomb[,"class"][[2]])
  # ref_class2 <- as.character(mat_gene_relcomb[,"class"][[1]])
  #vm(line, verbose)
  #print(ref_class1)
  #print(ref_class2)
  bam <- check_for_overlapping_reads(bamDna,
                                     bamRna,
                                     ref_chr1,
                                     ref_chr2,
                                     ref_pos1,
                                     ref_pos2, 
                                     verbose)
  ##prind(1bam)
  ##prind(1ref_ref1)
  ##prind(1ref_ref2)
  if(!length(bam)==0){  
    #vm("reads detected, applying core function", verbose, 1)
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
                               ref_class2,
                               verbose) %>%
      bind_rows() #%>%
    #  rowwise() %>%
    #  mutate(baseq1_conv=ascii_to_dec(baseq1),
    #         baseq2_conv=ascii_to_dec(baseq2))
    
    #vm("reads classified", verbose, 1)
    
  } else {
    #vm("no reads detected, return empty df", verbose, 1)
    classified_reads <- tibble()
  } 
  func_end(verbose)
  return(classified_reads)
}
phase_combination <- function(mat_gene_relcomb, comb, bamDna, bamRna, verbose, geneDir, 
                              phasing_type, showReadDetail){
  func_start()
  #print(mat_gene_relcomb)
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
                                          bamDna, bamRna, verbose)
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
                                                 printLog,
                                                 verbose
    ) %>%  
      mutate(
        dist=mat_dist[comb],
        class_comb=paste(mat_gene_relcomb[,"class"][[2]], mat_gene_relcomb[,"class"][[1]], 
                         sep="-"),
        comb=comb,
        phasing=comb
      ) 
  }   
  append_matrices(classified_main_comb)
  func_end()
  return(classified_main_comb)
}
append_matrices <- function(classified_main_comb, iterate=TRUE){
  func_start()
  mat_old <- mat_phased
  lapply(seq(1,nrow(classified_main_comb)), function(C){
    #print(classified_main_comb[C,]$comb)
    #print(classified_main_comb[C,]$nstatus)
    #split_list <- str_split(classified_main_comb[C,]$comb, "-") %>% unlist()
    #print(split_list)
    #mat_phased[split_list[1], split_list[2]] <<- classified_main_comb[C,]$nstatus
    #mat_conf[split_list[1], split_list[2]] <<- classified_main_comb[C,]$conf
    mat_phased[classified_main_comb[C,]$comb] <<- classified_main_comb[C,]$nstatus
    ##print(mat_phased)
    ##print(length(mat_phased))
    ##print(dim(mat_phased))
    ##print(nrow(mat_phased))
    ##print(ncol(mat_phased))
    #mat_phased[classified_main_comb[C,]$comb] <<- 1
    #mat_conf[classified_main_comb[C,]$comb] <<- classified_main_comb[C,]$conf
    mat_info[classified_main_comb[C,]$comb] <<- classified_main_comb[C,]$phasing
    #print(mat_phased)
    return()
  })
  if(iterate==TRUE){
    something_changed <- sum(mat_old, na.rm = T)!=sum(mat_phased, na.rm = T)
    mat_new <- mat_phased
    mat_new[is.na(mat_new)] <- 0
    mat_old[is.na(mat_old)] <- 0
    changes <- which(!mat_new==mat_old) 
    i <- 1
    ##solved_master_combs <<- list()
    while(i<=length(changes)){
     #message("left changes:", changes)
      ## something changed
      ##print("something changed")
      ## get mut/snp ids that have changed
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
         #print("state 1")
          conns <- sort(unique(c(m1_same_conn, m2_same_conn)))
          for(x in 1:(length(conns)-1)){
            #print(i)
            for(j in (x+1):length(conns)){
              ##print(j)
             #message(paste(conns[x],conns[j], sep="-"))
              p1 <- conns[x]
              p2 <- conns[j]
             #print(p1)
             #print(p2)
             #print(m12)
             #print(which(!m12 %in% c(p1,p2)))
              
              notxf <- m12[which(!m12 %in% c(p1,p2))]
              #print(notx)
             #print(111111)
             #print(notxf)
             #print(length(notxf))
              
              
              
              if(length(notxf)>=1){
                notx <- min(notxf)
               ##print(mat_phased[p1,p2])
                if(is.na(mat_phased[p1,p2])|mat_phased[p1,p2]==0){
                 #print("change mat_info")
                  #print(mat_info)
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
                  #print(combined_info)
                }
              } 
              if(is.na(mat_phased[p1,p2])|mat_phased[p1,p2]==0){
                mat_phased[p1,p2] <<- 1
              }
            }
          }
        } else {
         #print("state 2")
          ## if status == diff
          for(x in m1_same_conn){
            #print(i)
            for(j in m2_same_conn){
              ##print(j)
              #print(paste(conns[i],conns[j], sep="-"))
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
  func_end()
  #print(mat_phased)
}

perform_direct_phasing <- function(mat_gene, bamDna, bamRna, purity, 
                                   verbose, printLog, showReadDetail, geneDir){
  func_start()
  #phasingDir <- NULL
  phasing_type <- "direct"
  # if(!is.null(geneDir)){
  #   phasingDir <- file.path(geneDir, phasing_type)
  #   dir.create(phasingDir)
  # }
 #print(mat_phased)
 #print(mat_info)
  phasing_info <- tibble()
  unphased <- which(is.na(mat_phased))
  while(length(unphased)>0){
    #print(unphased)
    ## as long as there are unphased combinations, try to phase
    relcombxy <- get_xy_index(unphased[1], nrow(mat_phased))
    mat_gene_relcomb <- mat_gene[relcombxy,]
    comb <- unphased[1]
    classified_main_comb <- phase_combination(mat_gene_relcomb, comb, bamDna, bamRna, 
                                              verbose, geneDir, phasing_type, showReadDetail)
   #print(mat_phased)
   #print(mat_info)
    phasing_info <- bind_rows(phasing_info,
                              classified_main_comb)
    unphased <- which(is.na(mat_phased))
  }

  func_end()
  return(list(status=NULL,
              info=phasing_info,
              exit=tibble(phasing="direct", info="test")))
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
phase <- function(df_gene, bamDna, bamRna, 
                  showReadDetail, purity, sex, haploBlocks, phasedVcf,
                  distCutOff, printLog, verbose, logDir, refGen, somCna, 
                  snpQualityCutOff, 
                  phasingMode){
  vm(as.character(sys.call()[1]), verbose, 1)
  #haploblock_phasing <- NULL
  GENE <- unique(df_gene$gene)
  if(!is.null(logDir)){
    geneDir <- file.path(logDir, GENE)
    dir.create(geneDir)
  } else {
    geneDir <- NULL
  }
  ##prind(1df_gene)
  ## (1): define all combinations of variants to be phased
  #all_combinations <- make_phasing_combinations(df_gene) 
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
                                           purity, verbose, printLog, 
                                           showReadDetail, geneDir)
  phasing_info <- direct_phasing$info
  #print(phasedVcf)
  if(length(unknown_main())>0&!is.null(phasedVcf)){
    append_loglist("unphased combinations left --> Initialize SNP phasing")
    ## missing combinations in main muts --> start secondary phasing approaches
    lsnps <- load_snps(df_gene, phasedVcf, haploBlocks, refGen, distCutOff, verbose, somCna, snpQualityCutOff)
    store_log(geneDir, lsnps, "lsnps.tsv")
    #print(lsnps)
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
     #print(mat_phased)
      append_phasing_matrices(snps$mut_id, snps$start, distCutOff)
     #print(mat_phased)
     #print(snps)
      add_snps_to_matrices(snps)
     #print(mat_phased)
      ## build main daatframe of muts and snps and annotate haploblocks and segments of allelic imbalance
      df <- bind_rows(
        df_gene %>% select(chr, pos, ref, alt, class, af, mut_id),
        snps %>% select(chr=seqnames, pos=start, ref=REF, alt=ALT, af, mut_id, gt=gt_final) %>% mutate(class="snp")
      )


      to_phase <- prioritize_combination()
      i <- 1
      while(!is.null(to_phase)&length(unknown_main())>0){
        #print(to_phase)
        #print(unknown_main())
        append_loglist("Phasing combination:", paste(paste0(to_phase), collapse="-"))

        mat_gene_relcomb <- df[to_phase,] %>% as.matrix()
        comb <- get_single_index(to_phase[1],
                                 to_phase[2],
                                 nrow(mat_phased))

        classified_main_comb <- phase_combination(mat_gene_relcomb, comb,
                                                  bamDna, bamRna, verbose, 
                                                  geneDir, comb, showReadDetail)
        #print(classified_main_comb)
        phasing_info <- bind_rows(phasing_info,
                                  classified_main_comb)
        #print(main_muts)
        #print(snps$mut_id)
        to_phase <- prioritize_combination()
        i <- i +1
        #print(to_phase)
        #print(mat_phased[c(1:4), c(1:4)])
      }
      #print(df)
      store_log(geneDir, df, "df_snps_muts.tsv")
      #print(phasing_info)
      store_log(geneDir, phasing_info, "all_indirect_phasing_combinations.tsv")
    }## no snps found
  }
  ## secondary phasing done
  append_loglist("finalizing phasing results")


  ## reconstruct phasing results
  main_pos <- get_main_mut_pos()
  uppertri <- which(upper.tri(mat_phased))
  #solved <- main_pos[which(mat_phased[main_pos]>0)]
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
                  c(colnames(mat_phased)[get_xy_index(ncomb, nrow(mat_phased))["x"]], 
                     rownames(mat_phased)[get_xy_index(ncomb, nrow(mat_phased))["y"]]),
                  levels=c(paste0("m", seq(1,100,1)),
                           paste0("s", seq(1,2000,1)))
                 )
               ),
               collapse="-"
               ))   
 }

  
  store_log(geneDir, as_tibble(mat_phased), "mat_phased.tsv")
  #print(3)
  #print(as_tibble(mat_info))
  store_log(geneDir, as_tibble(mat_info), "mat_info.tsv")
  #print(4)
  store_log(geneDir, phasing_info, "all_phasing_combinations.tsv")
  func_end()
  return(list(phasing_all_combs, phasing_info_export, mat_phased, mat_info))
  
}