#' @keywords internal
catt <- function(printLog, level, text){
  if(printLog==TRUE){
    message(rep("  ", level), text)
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom GenomicAlignments readGAlignmentPairs first last
#' @importFrom Rsamtools ScanBamParam
#' @importFrom IRanges subsetByOverlaps
#' @importFrom dplyr tibble
prepare_raw_bam_file <- function(bamDna, chr1, chr2, pos1, pos2){
  qname.first <- . <- NULL
  ## importFrom dplyr tibble filter
  ref_pos1 <- as.numeric(pos1)
  ref_pos2 <- as.numeric(pos2)
  ref_chr1 <- as.character(chr1)
  ref_chr2 <- as.character(chr2)
  if(pos1>pos2){
    ref_pos1 <- as.numeric(pos2)
    ref_pos2 <- as.numeric(pos1)
    ref_chr1 <- as.character(chr2)
    ref_chr2 <- as.character(chr1)
  }
  ref_gr1 <- GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1)
  ref_gr2 <- GRanges(seqnames = ref_chr2, 
                                   ranges = ref_pos2)
  ## now load all reads/read-pairs that cover the position of the first variant
  all_covering_read_pairs <- readGAlignmentPairs(
    bamDna,
    param=ScanBamParam(
      which=GRanges(seqnames = ref_chr1, 
                                   ranges = ref_pos1),
      what=c("qname","seq", "cigar")
    )) 
  if(length(all_covering_read_pairs)==0){
    return(tibble())
  } else {
    ## combine all ranges and check for ref_pos2
    all_reads <- c(
      first(all_covering_read_pairs) %>%
        GRanges(),
      last(all_covering_read_pairs) %>%
        GRanges()
    )
    
    shared_read_pairs <- all_reads %>%
      subsetByOverlaps(.,ref_gr2) %>%
      elementMetadata(.) %>%
      .[["qname"]] %>%
      unique()
    
    filtered_reads <- all_reads[
      which(all_reads$qname %in% shared_read_pairs)
    ]
    return(filtered_reads)
  }
}
#' @keywords internal
allowed_inputs <- function(which_one){
  . <- NULL
  type_list <- list(
    cna_homdel_annotation = paste(c("HOMDEL", "HomoDel", "HomoDEL", "HOMODEL", 
                                    "HomDel", "homdel"),collapse = "|"),
    cna_incompletedel_annotation = c("DEL", "del", "Del"),
    chrom_names = c(
      c(seq_len(23), "X", "Y"),
      paste0("chr", c(seq_len(23), "X", "Y"))
    )
  )
  return(type_list[[which_one]])
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom purrr compact
#' @importFrom dplyr summarize pull
insert_missing_cnv_regions <- function(somCna, geneModel, sex, ploidy){
  seqnames <- start <- end <- tcn <- . <- NULL
  if(sex=="male"){
    allowed_chr <- allowed_inputs("chrom_names")
  } else {
    allowed_chr <- allowed_inputs("chrom_names") %>% 
      .[which(!str_detect(., "Y"))]
  }
  tbl_refgen <-
    as_tibble(geneModel) %>%
    group_by(seqnames) %>%
    summarize(mn=min(start),
              mx=max(end)) %>%
    filter(seqnames %in% allowed_chr)
  tbl_cnv <- 
    as_tibble(somCna) %>%
    mutate(tcn=as.numeric(tcn))
  new_somCna <- lapply(as.character(tbl_refgen$seqnames), function(CHR){
    all_chr <- tbl_cnv %>%
      filter(seqnames==CHR) %>%
      arrange(start)
    mn <- tbl_refgen %>% 
      filter(seqnames==CHR) %>% pull(mn)
    mx <- tbl_refgen %>% 
      filter(seqnames==CHR) %>% pull(mx)    
    
    if(nrow(all_chr)==0){
      tbl_refgen %>% 
        filter(seqnames==CHR) %>%
        select(seqnames, start=mn, end=mx) %>% 
        mutate(tcn=ploidy, cna_type="no CNA type annotated") %>%
        GRanges() %>%
        return()
    } else {
      if(all_chr[1,]$start<mn){mn <- all_chr[1,]$start}
      if(all_chr[nrow(all_chr),]$end>mx){mx <- all_chr[nrow(all_chr),]$end}
      poslist <- 
        apply(all_chr, 1, function(ROW){
          return(c(start=as.numeric(ROW[["start"]]), 
                   end=as.numeric(ROW[["end"]])))
        }, simplify=FALSE) %>%
        append(list(c(start=NA, end=mn)),.) %>%
        append(.,list(c(start=mx, end=NA)))
      
      to_insert <- lapply(c(2:length(poslist)), function(I){
        diff <- poslist[[I]][["start"]]-poslist[[I-1]][["end"]]
        if(diff>1){
          return(c(start=poslist[[I-1]][["end"]]+1, 
                   end=poslist[[I-1]][["end"]]+diff-1))
        }
      })  %>%
        compact()
      
      if(length(to_insert)==0){
        GRanges(all_chr) %>% return()
      } else {
        return(c(GRanges(all_chr),  bind_rows(to_insert) %>%
                   mutate(chr=CHR, tcn=ploidy, 
                          cna_type="no CNA type annotated") %>%
                   GRanges()))        
      }
    }
  })
  comb_somCna <- Reduce(function(x,y)c(x,y),new_somCna)
  return(comb_somCna)
}
pre_scoring <- function(tcn1, tcn2, status, aff_copies1, aff_copies2){
  mtcn <- min(as.numeric(tcn1), as.numeric(tcn2))
  left_wt_copies <- ifelse(status=="diff",
                           mtcn-sum(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)),
                           mtcn-max(as.numeric(aff_copies1), 
                                    as.numeric(aff_copies2)))
  pre_score <- ifelse(status=="diff"&left_wt_copies<0.5,
                      2, 1)
  return(list(mtcn=mtcn,
              left_wt_copies=left_wt_copies,
              pre_score=pre_score))
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr as_tibble group_by mutate tibble tally 
#' @importFrom purrr set_names
classify_combination <- function(classified_reads, purity, eval_full, printLog,
                                 tcn1, tcn2, aff_copies1, 
                                 aff_copies2){
  # tcn1 <- main_comb[['tcn1']]
  # tcn1 <- main_comb[['tcn2']]
  # aff_copies1<- main_comb[['aff_cp1']]
  # aff_copies2 <-  main_comb[['aff_cp2']]
  # classified_reads <- main_classified_reads
  result <- . <- fac <- NULL
  if(nrow(classified_reads)==0){
    return(rep(0, 7) %>% t() %>% as.data.frame() %>% 
             set_names(c("both", "none", "mut1", "mut2", "dev_var", 
                         "no_overlap", "none_raw")) %>% 
             mutate(status="null"))
  }
  number <- classified_reads %>%
    mutate(fac=factor(result, levels = c('both', 'mut1', 'mut2', 'none', 
                                         "dev_var", 'read_in_read', 
                                         'spanned_out'))) %>%
    group_by(fac, .drop = FALSE) %>%
    tally() %>% 
    column_to_rownames(var='fac') %>%
    t() %>%
    as_tibble()
  both <- number[['both']]
  mut1 <- number[['mut1']]
  mut2 <- number[['mut2']]
  none_raw <- number[['none']]
  # now the purity is calculated out... actually all "none" detetions  count 
  # for mono-allelic inactivation.... but to get a more reliable result it 
  # would be necessary to detect reads which have both mutations
  # therefore we check if any "both" detections were found and if not, we 
  # assume that the "none" detections are artefacts or non tumor cell tissue   
  none_adj_im <- none_raw-sum(number[1,])*(1-as.numeric(purity))
  none_adj <- ifelse(none_adj_im<0,
                     0,
                     round(none_adj_im))
  clear_sum <- sum(none_adj, both, mut1, mut2)
  prob_diff <- sum(mut1, mut2)/clear_sum
  prob_same <- sum(none_adj, both)/clear_sum
  ## here the main decision is made
  if(clear_sum==0){
    status <- "null"
    status_def_info <- "non-overlapping reads only -> null"
  } else if(sum(mut1, mut2, both)==0){
    ## if no evidence for both muts at the same read or at different reads can
    ## be found, we assign status "null"
    status <- "null"
    status_def_info <- 
      "Only reads without any variant: No evidence for same or diff -> null"
  } else if(eval_full==FALSE&none_raw==0&(mut1==0|mut2==0)){
    ## if we are at extended SNP phasing, indicated by the eval_full parameter
    ## we have to check if the variant we use to phase is a homozygous germline
    ## variant. In this case, we cannot use it for phasing, as it will be 
    ## detected at every read. To filter out such variants, we check if there 
    ## are any reads without the variant. If we cannot find any
    ## read in SNP phasing not having the variant, we assign "null" status as
    ## we assume the variant to phase is homozygous germline. this method is
    ## conservative.
    status <- "null"
    status_def_info <- paste(
      "one SNP of combination seems to be germline homozygous:",
      "Not usable for phasing -> null")
  } else if(both==0){
    status <- "diff"
    status_def_info <- ifelse(mut1==0|mut2==0,
      paste("no reads with both variants and reads with only one -> diff;",
      "but one variant without any detection!"),
      "variants found only on different reads -> diff"
    )
  } else if(prob_same>prob_diff){
    status <- "same"
    status_def_info <- ifelse(prob_diff>0,
      paste("Reads detected carrying both variants -> same;",
            "but also reads with only one!", "prob same:", prob_same, 
            "prob diff:", prob_diff),
      "variants detected on the same reads -> same"
      )
  } else if((mut1<0.1*sum(mut1, mut2)|mut2<0.1*sum(mut1, mut2))&
            both>0.1*sum(mut1, mut2)){
    status <- "same"
    status_def_info <- 
      paste("reads found with both variants and also with only one",
            "number of detections of one variant smaller than 10 % of",
            "overall evidence for diff and evidence for same higher -> same")
  } else {
    status <- "diff"
    status_def_info <- "variants found only on different reads -> diff"
  }
  ntabs <- ifelse(eval_full==TRUE, 3, 7)
  catt(printLog, ntabs, status_def_info)
  status_table <- tibble(
    both=both,
    none=none_adj,
    mut1=mut1,  
    mut2=mut2,
    dev_var=number[['dev_var']],
    no_overlap=sum(as.numeric(number[['read_in_read']]), 
                   as.numeric(number[['spanned_out']])),
    status=status,
    none_raw=none_raw
    )
  if(eval_full==TRUE){
    psc <- pre_scoring(tcn1, tcn2, status, aff_copies1, aff_copies2)
    if(status=="diff"&psc$pre_score==1&round(psc$mtcn)<=2){
      catt(printLog, 3, 
           "muts at different copies and tcn<=2, but left wt-copies >= 0.5!")
    }
    status_table <- status_table %>%
      mutate(
        score=psc$pre_score,
        wt_cp=round(psc$left_wt_copies,2),
        tcn=round(psc$mtcn,2),
        info=paste(status_def_info, "-> left wt copies", 
                   round(psc$left_wt_copies,2))
      )
  } 
  return(status_table)
}  
#' @keywords internal
#' @importFrom stringr %>% str_detect str_match_all str_sub
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate left_join bind_rows mutate tibble filter
make_read_overview <- function(read_start, seq, cigar){
  . <- element <- width <- end <- NULL
  raw_cigs <- str_match_all(cigar, 
                            '\\d*[:upper:]') %>% 
    unlist() %>%
    tibble(element=.) %>%
    mutate(
      width=str_sub(element, start = 1, end=nchar(element)-1) %>% 
        as.numeric(),
      type=str_sub(element, start = nchar(element), end=nchar(element))
    ) %>%
    filter(width!=0) %>%
    rownames_to_column("id")
  raw_maps <- 
    raw_cigs[which(!str_detect(raw_cigs$element, "D|N")),] %>%
    mutate(to_extract=cumsum(width)) %>%
    apply(.,1,function(cig){
      if(cig[["id"]]==1){
        curr_seq <- str_sub(seq, start=1, end=cig[["width"]])
      } else {
        curr_seq <- str_sub(seq, 
                            start=as.numeric(cig[["to_extract"]])-
                              as.numeric(cig[["width"]])+1, 
                            end=cig[["to_extract"]])
      }
      return(list(id=cig[["id"]],
                  seq=curr_seq))
    }) %>% bind_rows()
  comb <- left_join(raw_cigs, raw_maps, by="id")  
  rel <- 
    comb[which(comb$type %in% c("M", "D", "N")),] %>%
    mutate(end=cumsum(width),
           start=end-width+1,
           map_start=read_start+end-width,
           map_end=read_start+end-1)
  full <- left_join(comb, rel[c(1,8,9)], by="id")
  return(full)
}
#' @keywords internal
#' @importFrom stringr %>% str_replace str_sub
#' @importFrom dplyr between case_when rowwise filter
check_mut_presence <- function(read_structure, ref_pos,
                                  ref_alt, ref_ref, ref_class){
  map_start <- map_end <- . <- NULL
  ## check if read is spanned out:
  element_mut <- read_structure %>% rowwise() %>%
      filter(between(ref_pos, map_start, map_end)) %>%
      .[1,]
  if(element_mut$type %in% c("N", "D")){
    mut_at_read <- -10
  } else {
    if(ref_class=="snv"){
      pos_in_seq <- ref_pos-element_mut$map_start+1
      base <- 
        str_sub(element_mut$seq,
                start=pos_in_seq,
                end=pos_in_seq) 
      mut_at_read <- case_when(
        base==ref_alt ~ 1,
        base==ref_ref ~ 0,
        TRUE ~ -2
      )    
    } else if(ref_class=="ins"){
      cat_at_pos <- 
        read_structure[which(read_structure$map_end==ref_pos)+1,] %>%
        .[which(.$type=="I"),] 
      if(nrow(cat_at_pos)==0){
        mut_at_read <- 0
      } else {  
        length_ins <- str_replace(ref_alt, paste0("^", ref_ref), "") %>%
          nchar()
        if(length_ins==cat_at_pos[1,]$width){
          mut_at_read <- 1
        } else {
          mut_at_read <- -2
        }
      }
    } else if(ref_class=="del"){
      cat_at_pos <- 
        read_structure[which(read_structure$map_end==ref_pos)+1,]%>%
        .[which(.$type=="D"),]
      
      if(nrow(cat_at_pos)==0){
        mut_at_read <- 0
      } else {  
        length_ins <- str_replace(ref_ref, paste0("^", ref_alt), "") %>%
          nchar()
        if(length_ins==cat_at_pos[1,]$width){
          mut_at_read <- 1
        } else {
          mut_at_read <- -2
        }
      }
    }
  }
  return(mut_at_read)
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom purrr set_names
#' @importFrom dplyr case_when bind_rows as_tibble mutate
#' @importFrom GenomicRanges seqnames
core_tool <- function(qname, bam, 
                      ref_pos1, ref_pos2,
                      ref_alt1, ref_alt2, 
                      ref_ref1, ref_ref2,
                      ref_class1, ref_class2){
  . <- NULL
  paired_reads <- bam[which(bam$qname==qname)] %>%
    .[which(as.character(seqnames(.)) %in% 
              allowed_inputs("chrom_names"))] %>%
    as_tibble() %>%
    rownames_to_column("mate")
  if(length(paired_reads)!=0){
    
    read_structure <- apply(paired_reads, 1, function(READ){
      make_read_overview(as.numeric(READ[["start"]]),
                         READ[["seq"]],
                         READ[["cigar"]]) %>%
        mutate(mate=READ[["mate"]])
    }) %>% bind_rows() 
    mut1_in_read <- check_mut_presence(read_structure, ref_pos1,
                                          ref_alt1, ref_ref1, ref_class1)
    mut2_in_read <- check_mut_presence(read_structure, ref_pos2,
                                          ref_alt2, ref_ref2, ref_class2)
    final_assignment <- case_when(
      sum(mut1_in_read, mut2_in_read)==2 ~ "both",
      sum(mut1_in_read, mut2_in_read)==0 ~ "none",
      sum(mut1_in_read, mut2_in_read)<(-5) ~ "spanned_out",
      mut1_in_read==1 ~ "mut1",
      mut2_in_read==1 ~ "mut2",
      TRUE ~ "dev_var"
    ) 
  } else {
    ## non valid chromosomes of both reads
    return(NULL)
  }
  return(c(qname=qname, result=final_assignment, 
           origin=unique(paired_reads$origin)))
}
#' @keywords internal
#' description follows
#' @importFrom dplyr between
check_af <- function(af){
  num_af <- as.numeric(af)
  if(is.na(num_af)){
    stop("input allele frequency (af) is not numeric")
  } else if(!between(num_af, 0, 1)){
    stop("input allele frequency (af) must be between 0 and 1")
  } else {
    return(num_af)
  }
}
#' @keywords internal
#' description follows
check_tcn <- function(tcn){
  num_tcn <- as.numeric(tcn)
  if(is.na(num_tcn)){
    stop("input total copynumber (tcn) is not numeric")
  } else {
    return(num_tcn)
  }
}
#' @keywords internal
check_chr <- function(chr){
  if(chr %in% allowed_inputs("chrom_names")){
    return(as.character(chr))
  } else {
    stop("input chromosome (chr) must be either a number between ",
          "1 and 23 or X or Y, or in the chr1 format")
  }
}
#' calculates how many copies are affected by a germnline small variant
#'
#' @param af Allele-frequency of the variant (numeric value between 0 and 1)
#' @param tcn total-copynumber at position of the variant (numeric value >0)
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param chr chromosome of the variant (either format 1,2,..,X,Y or 
#' chr1,..,chrX)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param c_normal expected copy number at position of the variant in normal 
#' tissue, 1 for gonosomes in male samples, and 2 for male autosomes and all 
#' chromosomes in female samples. (The function can also assess the c_normal 
#' parameter by itself, but then the following two inputs must be provided: 
#' chr and sex)
#' @return A numeric value indicating the affecting copies for the variant
#' @export
#' @examples 
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#' aff_germ_copies(af=0.67, tcn=2, purity=0.9, chr="chrX", sex="female")
aff_germ_copies <- function(chr, af, tcn, purity, sex, 
                            c_normal=NULL){
  cv <- formula_checks(chr, af, tcn, purity, sex, c_normal)
  aff_cp <- cv$af*(cv$tcn+cv$c_normal*(1/cv$purity - 1)) - 1/cv$purity + 1
  return(aff_cp) 
  # alternative: (af*(purity*tcn+c_normal*(1-purity))-cc*(1-purity))/purity 
}
#' calculates how many copies are affected by a somatic small variant
#'
#' @param af Allele-frequency of the variant (numeric value between 0 and 1)
#' @param tcn total-copynumber at position of the variant (numeric value >0)
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param chr chromosome of the variant (either format 1,2,..,X,Y or 
#' chr1,..,chrX)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param c_normal expected copy number at the position of the variant in 
#' normal tissue, 1 for gonosomes in male samples, and 2 for male autosomes and 
#' all chromosomes in female samples. (The function can also assess the 
#' c_normal parameter by itself, but then the following two inputs must be 
#' provided: chr and sex)
#' @return A numeric value indicating the affecting copies for the variant
#' @examples 
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#' aff_som_copies(chr="chrX", af=0.67, tcn=2, purity=0.9, sex="female")
#' @export
aff_som_copies <- function(chr, af, tcn, purity, sex, c_normal=NULL){
  cv <- formula_checks(chr, af, tcn, purity, sex, c_normal)
  aff_cp <- cv$af*(cv$tcn+cv$c_normal*((1/cv$purity)-1))
  return(aff_cp)
  ## alternative: aff_cp = af*(tcn+c_normal/purity-c_normal)
}
#' @keywords internal
#' description follows
formula_checks <- function(chr, af, tcn, purity, sex, c_normal){
  purity <- check_purity(purity)
  af <- check_af(af)
  tcn <- check_tcn(tcn)
  if(is.null(c_normal)){
    sex <- check_sex(sex)
    chr <- check_chr(chr)
    if((sex=="male"&(chr=="X"|chr=="chrX"))|(chr=="Y"|chr=="chrY")){
      c_normal <- 1
    } else {
      c_normal <- 2
    }
  } else {
    c_normal <- check_ploidy(c_normal)
  }
  return(list(af=af, tcn=tcn, purity=purity, c_normal=c_normal))
}
#' @keywords internal
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate select mutate_at rowwise
#' @importFrom stringr str_detect
merge_sCNAs <- function(obj, somCna){
  alt <- af <- tcn <- tcn_assumed <- cna_type <- gene <- ref <- NULL
  obj %>%
    mergeByOverlaps(somCna) %>% as_tibble() %>%
    mutate(cna_type=ifelse(str_detect(cna_type, "LOH"), "LOH", "HZ")) %>%
    select(chr=1, pos=2, gene, ref, alt, af, tcn, cna_type, tcn_assumed) %>%
    mutate_at(.vars = c("af", "tcn"), .funs=as.numeric) %>%
    rowwise() %>%
    return()
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges elementMetadata
#' @importFrom IRanges mergeByOverlaps
#' @importFrom dplyr mutate filter
prepare_somatic_variant_table <- function(somSmallVars, templateGenes, 
                                          somCna, purity, sex){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- 
    aff_cp <- wt_cp <- . <- tcn_assumed <- NULL
  if(!is.null(somSmallVars)){
    rel <- somSmallVars%>%
      merge_sCNAs(., somCna) %>%
      
      ################################
      mutate(
        origin="somatic",
        class=define_class(ref, alt),
        aff_cp = aff_som_copies(chr, af, tcn, purity, sex),
        wt_cp = tcn-aff_cp,
        pre_info = 
          #ifelse(base::isTRUE(str_detect(cna_type,'LOH')), 
          ifelse(isTRUE(str_detect(cna_type,'LOH')),
            ifelse(wt_cp<=0.5,
              paste(
                  'somatic-variant -> LOH -> left wt copies',round(wt_cp,2),
                  '-> all copies affected'),
              paste(
                  'somatic-variant -> LOH -> ',
                  "left wt copies",round(wt_cp,2),"-> wt copies left")),
            ifelse(chr=="X"&sex=="male",
              ifelse(wt_cp<=0.5,
                paste(
                  "somatic-variant -> chrX & male -> ",
                  "left wt copies",round(wt_cp,2),"-> all copies affected"),
                paste(
                  'somatic-variant -> chrX & male -> left wt copies ',
                  round(wt_cp,2),
                  '-> wt copies left')),
                'somatic-variant -> no LOH -> wt copies left')
        )
      )%>%
      filter(gene %in% templateGenes)
    if(nrow(rel)==0){
      return(NULL)
    } else {
      return(rel)
    }
  } else {
    return(NULL)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr compact
#' @importFrom dplyr filter as_tibble select bind_rows mutate
extract_all_dels_of_sample <- function(somCna, geneModel, DEL_TYPE, 
                                       byTcn, sex, ploidy, include_dels){
  tcn <- cna_type <- seqnames <- gene <- . <- NULL
  if(include_dels==TRUE){
    TBL_CNV <- somCna %>%
      as_tibble() %>%
      mutate(tcn=as.numeric(tcn))
    if(DEL_TYPE=="homdel"){
      if(byTcn==TRUE){
        relevant_CNAs <- TBL_CNV %>%
          filter(round(tcn)==0)
      } else {
        relevant_CNAs <- TBL_CNV %>%
          filter(str_detect(cna_type, allowed_inputs("cna_homdel_annotation"))) 
      }
      pre_info <- 'homdel -> all copies affected'
    } else {
      if(byTcn==TRUE){
        if(sex=="female"){
          relevant_CNAs <- TBL_CNV %>%
            filter(round(tcn)<ploidy)
        } else {
          relevant_CNAs <- TBL_CNV %>%
            filter(!str_detect(seqnames, "X|Y")&round(tcn)<ploidy)
        }
      } else {
        relevant_CNAs <- TBL_CNV %>%
          filter(!str_detect(cna_type, 
                             allowed_inputs("cna_homdel_annotation"))) %>%
          filter(str_detect(cna_type, 
                            allowed_inputs("cna_incompletedel_annotation")))
      }
      pre_info <- "incompletedel -> wt-copies left"
    } 
    if(nrow(relevant_CNAs)==0){
      return(NULL)
    } else {
      merged_CNAs <- apply(relevant_CNAs, 1, function(CNV){
        GR_single_CNV <- GRanges(as_tibble(t(CNV)))
        suppressWarnings(
        rel <- subsetByOverlaps(geneModel, GR_single_CNV) %>%
          as_tibble()
        )
        if(nrow(rel)==0){
          return(NULL)
        } else {
          rel %>%
            select(gene) %>%
            mutate(chr=CNV[["seqnames"]],
                   pos=as.numeric(CNV[["start"]]),
                   tcn=as.numeric(CNV[["tcn"]]),
                   tcn_assumed=FALSE,
                   af=NA,
                   origin="somatic",
                   class=DEL_TYPE,
                   ref= NA,
                   alt= NA,
                   cna_type=DEL_TYPE,
                   aff_cp=ploidy-as.numeric(CNV[["tcn"]]),
                   wt_cp=as.numeric(CNV[["tcn"]]),
                   pre_info=pre_info,
            ) %>%
            return()
        }
      })%>% compact() %>% bind_rows()
      if(nrow(merged_CNAs)==0){
        return(NULL)
      } else {
        return(merged_CNAs)
      }
    }
  } else {
    return(NULL)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom dplyr mutate select
prepare_germline_variants <- function(germSmallVars, somCna, purity, sex){
  cna_type <- gene <- ref <- alt <- af <- tcn <- cna_type <- chr <- aff_cp <- 
    origin <- pos <- wt_cp <- pre_info <- . <- tcn_assumed <- NULL
  if(is.null(germSmallVars)){
    return(NULL)
  } else {
    df_germ <- germSmallVars %>% 
      merge_sCNAs(., somCna) %>%
      mutate(
        origin="germline",
        class=define_class(ref, alt),
        aff_cp = aff_germ_copies(chr, af, tcn, purity, sex, 1),
        wt_cp = tcn-aff_cp,
        pre_info = 
          ifelse(isTRUE(cna_type=='LOH'),
            ifelse(af>=0.5,
                 paste(
                   'germline-variant -> LOH ->',
                   'AF >= 0.5 -> all copies affected'),
                 paste(
                   'germline-variant -> LOH ->',
                   ' AF < 0.5 -> variant lost in tumor')),
            'germline-variant -> no LOH -> wt copies left')
      ) %>%
      select(gene,
             origin,
             class,
             chr,
             pos,
             ref,
             alt,
             af,
             tcn, 
             tcn_assumed,
             cna_type,
             aff_cp,
             wt_cp,
             pre_info
      ) 
    return(df_germ)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom GenomicRanges elementMetadata elementMetadata<-
#' @importFrom methods is
check_somCna <- function(somCna, geneModel, sex, ploidy,
                         assumeSomCnaGaps, colnameTcn, 
                         colnameCnaType){
  . <- NULL
  ## check if class i GRanges
  if(!is(somCna, "GRanges")){
    stop(
      "input somCna must be a GRanges object; given input appears to be: ", 
      class(somCna))
  } else if(
    !(
      ("tcn" %in% names(elementMetadata(somCna))|
       !is.null(colnameTcn))&
      ("cna_type" %in% names(elementMetadata(somCna))|
       !is.null(colnameCnaType))
    )
    ){
    stop("input somCna requires the following metadata columns: \'tcn\'",
             "and \'cna_type\'",
             "please rename relevant metadata or provide metadata colname by",
             "colnameTcn or colnameCnaType")
  } else {
    if(!is.null(colnameTcn)){
      elementMetadata(somCna)[,"tcn"] <- 
        elementMetadata(somCna)[,colnameTcn]
    }
    if(!is.null(colnameTcn)){
      elementMetadata(somCna)[,"cna_type"] <- 
        elementMetadata(somCna)[,colnameCnaType]
    }
    somCna$tcn_assumed <- FALSE
    if(sum(is.na(elementMetadata(somCna)[,"cna_type"]))>0){
      warning("cna_type column of input somCna contains", 
              sum(is.na(elementMetadata(somCna)[,"cna_type"])),
                    "NA values;\n  they will be taken as hetero-zygous\n")
      elementMetadata(somCna)[,"cna_type"][
        which(
          is.na(elementMetadata(somCna)[,"cna_type"]))] <- NA
    } 
    if(assumeSomCnaGaps==TRUE){
      if(sum(is.na(elementMetadata(somCna)[,"tcn"]))>0){
        warning("tcn column of input somCna contains", 
                      sum(
                        is.na(elementMetadata(somCna)[,"tcn"])),
                      "NA values;\n  they will be taken as ground ploidy:") 
        elementMetadata(somCna)[,"tcn_assumed"][
          which(is.na(elementMetadata(somCna)[,"tcn"]))] <- TRUE
        elementMetadata(somCna)[,"tcn"][
          which(is.na(elementMetadata(somCna)[,"tcn"]))] <- 
          ploidy
      }
      new_somCna <- insert_missing_cnv_regions(somCna, geneModel, 
                                               sex, ploidy)
    } else {
      if(sum(is.na(elementMetadata(somCna)[,"tcn"]))>0){
        warning("tcn column of input somCna contains", 
                      sum(
                        is.na(elementMetadata(somCna)[,"tcn"])),
                      "NA values;\n  they will be excluded from the analysis.",
                      "Use assumeSomCnaGaps=TRUE to take ground ploidy",
                      "(ploidy) as tcn"
                    )
        new_somCna <- somCna[which(!is.na(
          elementMetadata(somCna)[,"tcn"]))]
      } else {
        new_somCna <- somCna
      }
      
    }
    return(new_somCna)
  }
}
#' @keywords internal
#' description follows
#' @importFrom GenomicRanges elementMetadata
nm_md <- function(obj){
  return(names(elementMetadata(obj)))
}
#' @keywords internal
#' description follows
check_name_presence <- function(obj, type){
  if(
      (   
        ("gene" %in% nm_md(obj)|
         "GENE" %in% nm_md(obj)
        )&
        (type=="gene_model"|
          (
            ("af" %in% nm_md(obj)|
             "AF" %in% nm_md(obj)
            )&
            ("ref" %in% nm_md(obj)|
             "REF" %in% nm_md(obj)
            )&
            ("alt" %in% nm_md(obj)|
             "ALT" %in% nm_md(obj)
            )
          )
        )
      )
    ){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_match
#' @importFrom GenomicRanges elementMetadata elementMetadata<-
assign_correct_colnames <- function(obj, type){
  . <- NULL
  col_gene <- str_match(
    nm_md(obj), 
    "gene|GENE") %>%
    .[which(!is.na(.))]
  elementMetadata(obj)[,"gene"] <- 
    elementMetadata(obj)[,col_gene]  
  if(type=="small_vars"){
    col_af <- str_match(
      nm_md(obj), "af|AF") %>%
      .[which(!is.na(.))]
    col_ref <- str_match(
      nm_md(obj), "ref|REF") %>%
      .[which(!is.na(.))]
    col_alt <- str_match(
      nm_md(obj), "alt|ALT") %>%
      .[which(!is.na(.))]
    elementMetadata(obj)[,"af"] <- 
      elementMetadata(obj)[,col_af]
    elementMetadata(obj)[,"ref"] <- 
      elementMetadata(obj)[,col_ref]
    elementMetadata(obj)[,"alt"] <- 
      elementMetadata(obj)[,col_alt]
    elementMetadata(obj)["mid"] <- 
      seq_len(length(obj))    
  }
  return(obj)
}

general_gr_checks <- function(obj, type, lab){
  if(!is(obj, "GRanges")){
    stop("input ", lab, " must be a GRanges object;",
         "given input appears to be:", 
         class(obj))
  } else if(!check_name_presence(obj, type)){
    if(type=="small_vars"){
      stop("input ", lab, " requires the following metadata columns:",
           "\'gene\'/\'GENE\', \'ref\'/\'REF\',",
           " \'alt\'/\'ALT\' and \'af\'/\'AF\'")      
    } else {
      stop(
        "input geneModel requires the following metadata columns: ",
        "\'gene\'/\'GENE\'")
    }
  } else {
    return(assign_correct_colnames(obj, type))
  }
}

#' @keywords internal
#' description follows
check_gr_gene_model <- function(geneModel){
  . <- NULL
  if(is.null(geneModel)){
    stop("input geneModel must not be NULL")
  } else {
    return(general_gr_checks(geneModel, "gene_model", "geneModel"))
  }
}
#' @keywords internal
#' description follows
check_gr_small_vars <- function(obj, origin){
  . <- NULL
  lab <- ifelse(origin=="somatic",
                 "somSmallVars",
                 "germSmallVars")
  if(!is.null(obj)){
    return(general_gr_checks(obj, "small_vars", lab))
  } else {
    warning("Input ", lab, " empty/does not contain variants. ",
                  "Assuming there are no ", origin, " small variants")
    return(NULL)
  }
}
#' @keywords internal
#' description follows 
#' @importFrom dplyr between
check_purity <- function(purity){
  if(is.na(as.numeric(purity))){
    stop("input purity must be numeric or a character that can",
            "be converted to numeric;\n  ", purity, 
            "can not be converted to numeric")
  } else if(!between(as.numeric(purity), 0, 1)){
    stop("input purity must be a numeric value between 0 and 1")
  } else {
    return(as.numeric(purity))
  }
}
#' @keywords internal
#' description follows
check_ploidy <- function(ploidy){
  if(is.null(ploidy)){
    return(NULL)
  } else {
    if(is.na(as.numeric(ploidy))){
      stop("input ploidy/c_normal must be numeric or a character that can",
              "be converted to numeric;\n  ", ploidy, 
              "can not be converted to numeric")
    } else {
      return(as.numeric(ploidy))
    }
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_detect
check_sex <- function(sex){
  . <- NULL
  allowed_sex <- c("male", "m", "female", "f") %>% c(.,toupper(.))
  if(!sex %in% allowed_sex){
    allowed_sex_to_print <- paste(allowed_sex, collapse = "\', \'") %>% 
              paste0("\'", ., "\'")
    stop("input sex must be one of ", allowed_sex_to_print)
  } else {
    sex <- tolower(sex)
    if(str_detect(sex, "f")){
      return("female")
    } else {
      return("male")
    }
  }
}
#' @keywords internal
#' description follows 
check_bam <- function(bamDna){
  if(file.exists(bamDna)){
    return(bamDna)
  } else {
    stop("input bamDna does not exist")
  }
}
#' @keywords internal
#' description follows
check_vcf <- function(vcf){
  if(is.null(vcf)){
    return(NULL)
  } else if(sum(unlist(lapply(vcf, file.exists)))==length(vcf)){
    return(vcf)
  } else {
    stop("input vcf does not exist")
  }
}
#' @keywords internal
#' description follows
check_rna <- function(bamRna){
  if(is.null(bamRna)){
    message("no RNA file provided: Analysis will be done without RNA reads")
    return(NULL)
  } else if(file.exists(bamRna)){
    return(bamRna)   
  } else {
    warning("input bamRna does not exist:", bamRna,
                  "\nanalysis will be done without RNA reads\n")
    return(NULL)
  }
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_split
#' @importFrom purrr map_chr set_names
make_phasing_combinations <- function(df_mut_to_combine){
  . <- pos1 <- pos2 <- mut_id1 <- mut_id2 <- NULL
  outer(df_mut_to_combine$pos, df_mut_to_combine$pos, `paste`) %>%
    .[which(upper.tri(.))] %>%
    as.data.frame() %>% 
    set_names(nm='raw') %>%
    mutate(
      pos1=str_split(raw, ' ') %>%map_chr(.,1) %>% as.numeric(),
      pos2=str_split(raw, ' ') %>%map_chr(.,2) %>% as.numeric(),
      distance = abs(pos1-pos2)
    ) %>%
    # add all information about every mutation to the data frame
    merge(df_mut_to_combine %>% 
            set_names(nm=names(df_mut_to_combine) %>% paste0(.,'1')), 
          by='pos1', all.x=TRUE) %>%
    merge(df_mut_to_combine %>% 
            set_names(nm=names(df_mut_to_combine) %>% paste0(.,'2')), 
          by='pos2', all.x=TRUE) %>%
    .[order(.$distance),] %>%
    mutate(comb_id=paste(mut_id1, mut_id2, sep='-'))
}
#' @keywords internal
#' description follows
check_for_overlapping_reads <- function(bamDna, bamRna,
                                        ref_chr1, 
                                        ref_chr2, 
                                        ref_pos1, 
                                        ref_pos2){
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
  return(c(dna_bam, rna_bam))
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr tibble 
classify_reads <- function(line, bamDna, bamRna){
  dist <- abs(as.numeric(line[['pos1']])-as.numeric(line[['pos2']]))
  ref_pos1 <- as.numeric(line[['pos1']])
  ref_pos2 <- as.numeric(line[['pos2']])
  ref_chr1 <- as.character(line[['chr1']])
  ref_chr2 <- as.character(line[['chr2']])
  ref_alt1 <- as.character(line[['alt1']])
  ref_alt2 <- as.character(line[['alt2']])
  ref_ref1 <- as.character(line[['ref1']])
  ref_ref2 <- as.character(line[['ref2']])
  ref_class1 <- as.character(line[['class1']])
  ref_class2 <- as.character(line[['class2']])
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
  return(classified_reads)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom Rsamtools TabixFile
#' @importFrom VariantAnnotation readVcf
#' @importFrom DelayedArray rowRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom dplyr rename arrange as_tibble left_join mutate ungroup rowwise ungroup select
check_for_snps_between_main_muts <- function(main_comb, vcf, df_gene){
  . <- ALT <- REF <- end <- start <- seqnames <- chr <- pos <- mut_id <- NULL
  pos_minor <- min(main_comb[['pos2']],
                   main_comb[['pos1']])
  pos_major <- max(main_comb[['pos2']],
                   main_comb[['pos1']])
  gr_roi <- GRanges(paste0(main_comb[["chr1"]], ":",
                                          pos_minor, "-", pos_major))
  vars_in_between_raw <- lapply(vcf, function(VCF){
    tab_vcf <- TabixFile(VCF)
    vcf_region <- 
      readVcf(tab_vcf, "hg19", 
                                 param=gr_roi) %>%
      #DelayedArray::rowRanges() %>%
      rowRanges() %>%
      return()
  }) %>%
    Reduce(function(x,y)c(x,y),.) %>%
    as_tibble() %>%
    ## maybe improve handling of more entries in DNAstringLists
    ## now it takes the first
    mutate(ALT=unlist(lapply(ALT, function(x){as.character(x[[1]][1])}))) %>%
    rowwise() %>%
    mutate(ALT=as.character(ALT),
           class=define_class(REF, ALT)) %>%
    select(1,2,3, ref=REF, alt=ALT, class) %>%
    unique() %>%
    select(-end) %>%
    rename(pos=start, chr=seqnames) %>%
    #dplyr::rename(pos=start, chr=seqnames) %>%
    left_join(df_gene %>% select(chr, pos, mut_id),
              by=c("chr"="chr",
                   "pos"="pos")) %>%
    arrange(mut_id) %>%
    ungroup() %>%
    ## remove main mutations if they are in between... because would be
    ## duplicated with upper analysis
    filter(is.na(mut_id)|
             mut_id %in% c(main_comb[["mut_id1"]], main_comb[["mut_id2"]]))
  return(vars_in_between_raw)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr mutate as_tibble
check_read_presence <- function(sub_comb, bamDna, bamRna){
  bam_raw <- check_for_overlapping_reads(bamDna, bamRna,
                                         sub_comb[["chr1"]],
                                         sub_comb[["chr2"]],
                                         sub_comb[["pos1"]],
                                         sub_comb[["pos2"]]) 
  if(length(bam_raw)!=0){
    bam <- bam_raw 
    bam$comb_id <- sub_comb[["comb_id"]]
  } else {
    bam <- NULL
  }
  res_tibble <- as_tibble(t(sub_comb)) %>%
    mutate(nreads=0.5*length(bam_raw))
  rl <- list(rt=res_tibble, bam=bam)
  return(rl)
}
#' @keywords internal
find_path <- function(still_present_muts, sub_checked_read_presence,
                      main_comb){
  . <- mut_id1 <- mut_id2 <- NULL
  a_star_nodes <- lapply(still_present_muts, function(MUT){
    rel_rows <- sub_checked_read_presence %>%
      filter(mut_id1==MUT|mut_id2==MUT)
    all_partners <- c(rel_rows$mut_id1, rel_rows$mut_id2) %>%
      .[which(!.==MUT)]
    return(lapply(all_partners, function(x)return(1)) %>% 
             unlist())
  })
  path <- a_star_pathfinder(
    a_star_nodes,
    main_comb[["mut_id1"]],
    main_comb[["mut_id2"]])
  return(path)
}
#' @keywords internal
#' @importFrom stringr %>% str_split
#' @importFrom dplyr bind_rows mutate_at
#' @importFrom purrr compact
phase_along_path <- function(path, sub_checked_read_presence, bam_raw,
                             purity, printLog){
  . <- mut_id1 <- mut_id2 <- comb_id <- NULL
  rel_sub_combs <- lapply(seq_len(length(path)), function(N){
    if(N==length(path)){
      return(NULL)
    } else {
      sub_checked_read_presence %>%
        filter(
          (path[[N]]==mut_id1|path[[N]]==mut_id2)&
            (path[[N+1]]==mut_id1|path[[N+1]]==mut_id2)
        ) %>%
        return()
    }
  }) %>%
    compact() %>%
    bind_rows() %>%
    mutate_at(.vars = c("chr1", "chr2", "pos1", "pos2", "distance"),
              .funs=as.numeric)
  bam <- bam_raw %>% compact() %>% Reduce(function(x,y)c(x,y),.)
  sub_fin_comb <-  apply(rel_sub_combs, 1, function(sub_comb){
    catt(printLog, 6, c("sub combination:", sub_comb[["comb_id"]], 
                        "distance:", sub_comb[["distance"]]))
    bam_fil <- bam %>% .[which(.$comb_id == sub_comb[["comb_id"]])]
    sub_classified_reads <- lapply(unique(bam_fil$qname), 
                                  core_tool,
                                  bam_fil,
                                  as.numeric(sub_comb[["pos1"]]),
                                  as.numeric(sub_comb[["pos2"]]),
                                  sub_comb[["alt1"]],
                                  sub_comb[["alt2"]],
                                  sub_comb[["ref1"]],
                                  sub_comb[["ref2"]],
                                  sub_comb[["class1"]],
                                  sub_comb[["class2"]]) %>%
      compact() %>%
      bind_rows()
    #print(sub_classified_reads)
    sub_status_combination <- 
      classify_combination(sub_classified_reads,
                           purity,
                           eval_full=FALSE,
                           printLog)
    return(sub_status_combination %>%
             mutate(comb_id=sub_comb[["comb_id"]]))
  }) %>% bind_rows()%>%
    mutate(
      f=str_split(comb_id, "-") %>% map_chr(.,1),
      l=str_split(comb_id, "-") %>% map_chr(.,2),
    )
  return(sub_fin_comb)
}
#' @keywords internal
#' @importFrom stringr %>% str_count str_replace_all str_detect
#' @importFrom dplyr desc
recombine_phasing_results <- function(non_null, main_comb){
  status <- comb_id <- f <- l <- fin_comb <- NULL
  sames <- 
    non_null %>%
    arrange(desc(status))  %>%
    select(status, comb_id, f, l)
  new <- sames
  for(n in seq_len(str_count(sames$status, "same") %>% sum())){
    curr_comb_id <- new[[n, "comb_id"]]
    to_rep <- paste(paste0(new[[n, "f"]],"$"),
                    paste0(new[[n, "l"]],"$"),
                    sep="|")
    replacement <- paste0(new[[n, "f"]],new[[n, "l"]])
    new <- new %>% 
      mutate(f=str_replace_all(f, to_rep, replacement),
             l=str_replace_all(l, to_rep, replacement)) 
  }
  final <- new %>%
    mutate(fin_comb=paste(f, l, sep="-")) %>%
    filter(str_detect(fin_comb, main_comb[["mut_id1"]])&
             str_detect(fin_comb, main_comb[["mut_id2"]]))
  return(final)
}
#' @keywords internal
#' @importFrom dplyr tibble
finalize_snp_phasing <- function(final, main_comb){
  status <- unique(final$status)
  tcn1 <- main_comb[["tcn1"]]
  tcn2 <- main_comb[["tcn2"]]
  aff_copies1 <- main_comb[["aff_cp1"]]
  aff_copies2 <- main_comb[["aff_cp2"]]
  psc <- pre_scoring(tcn1, tcn2, status, aff_copies1, aff_copies2)
  RESULT <- tibble(
    comb=main_comb[['comb_id']],
    class_comb=paste(main_comb[['class1']], 
                     main_comb[['class2']], sep="-"),
    dist=as.numeric(main_comb[["distance"]]),
    score=psc$pre_score,
    wt_cp=round(psc$left_wt_copies,2),
    tcn=round(psc$mtcn,2),
    info="status defined after extended HP phasing",
    status=status,
    DNA_rds=0,
    RNA_rds=0,
    both=0,
    none=0,
    mut1=0,
    mut2=0,
    no_overlap=0,
    none_raw=0,
    dev_var=0
  )
  return(RESULT)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr mutate
return_null_result <- function(main_comb, bamRna, ext_snp_info){
  RESULT <- rep(0, 9) %>% t() %>% as.data.frame() %>% 
    set_names(c("both", "none", "mut1", "mut2", "dev_var", 
                "no_overlap", "none_raw", "DNA_rds", "RNA_rds")) %>% 
    mutate(
      comb=main_comb[['comb_id']],
      class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                       sep="-"),
      dist=as.numeric(main_comb[["distance"]]),
      status="null",
      score=1,
      wt_cp=min(as.numeric(main_comb[['tcn1']]), 
                as.numeric(main_comb[['tcn2']])) %>% 
        round(2),
      tcn=min(as.numeric(main_comb[['tcn1']]), 
              as.numeric(main_comb[['tcn2']])) %>% 
        round(2),
      info=ifelse(!is.null(bamRna),
                  paste('unclear: no ovrlp. reads;', ext_snp_info),
                  paste('unclear: no ovrlp. reads & no RNA file;', 
                        ext_snp_info))
    )
} 
#' @keywords internal
#' description follows
#' @importFrom stringr %>%
#' @importFrom purrr map_chr set_names
#' @importFrom dplyr mutate nth select filter bind_rows
phase <- function(df_gene, bamDna, bamRna, 
                  showReadDetail, purity, vcf, distCutOff, printLog){
  pos1 <- pos2 <- mut_id1 <- mut_id2 <- qname.first <- result <- . <- origin <- 
    qname <- mut_id <- nreads <- status <- NULL
  # now a dataframe with all necesarry combination of input mutations is created
  all_combinations <- make_phasing_combinations(df_gene) 
  #main_comb <- all_combinations[1,]
  final_combinations <-  apply(all_combinations, 1, function(main_comb){
    catt(printLog, 2 ,c("combination:", main_comb[["comb_id"]]))
    return_null_result <- TRUE
    try_snp_phasing <- TRUE
    ext_snp_phasing <- NULL
    main_classified_reads <- classify_reads(main_comb, bamDna, bamRna)
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
                                     eval_full = TRUE, 
                                     printLog,
                                     main_comb[['tcn1']],
                                     main_comb[['tcn2']],
                                     main_comb[['aff_cp1']],
                                     main_comb[['aff_cp2']]
                                     ) %>%  
        mutate(
          DNA_rds=nrow(main_classified_reads %>% filter(origin=="DNA")),
          RNA_rds=nrow(main_classified_reads %>% filter(origin=="RNA")),
          dist=as.numeric(main_comb[["distance"]]),
          class_comb=paste(main_comb[['class1']], main_comb[['class2']], 
                           sep="-"),
          comb=main_comb[['comb_id']]
        ) 
      if(RESULT$status!="null"){
        try_snp_phasing <- FALSE
        return_null_result <- FALSE
        catt(printLog, 3, c("status:", RESULT$status))
      } else {
        catt(printLog, 3, "status can not be defined")
      }
    }   
    if(is.null(vcf)){
      ext_snp_info <- "no vcf provided for extended SNP phasing"
    } else if(try_snp_phasing==FALSE){
      ext_snp_info <- ""
    } else if(as.numeric(main_comb[["distance"]])<distCutOff){
      catt(printLog, 4, "Initializing extended SNP phasing")
      vars_in_between_raw <- check_for_snps_between_main_muts(main_comb, vcf,
                                                              df_gene)
      if(nrow(vars_in_between_raw)>2){
        catt(printLog, 5, c(nrow(vars_in_between_raw)-2, 
              "SNPs between variants detected"))
        vars_in_between <- vars_in_between_raw %>% 
          mutate(mut_id=
                   c(mut_id[which(!is.na(mut_id))],
                     paste0("s", 
                            seq_len(nrow(.)-length(which(!is.na(mut_id)))))))
        sub_combinations <- make_phasing_combinations(vars_in_between) %>%
          .[which(!(str_detect(.$mut_id1, "m")&str_detect(.$mut_id2, "m"))),]
        sub_checked_read_presence_list <- apply(sub_combinations, 1, 
                                                check_read_presence, 
                                                bamDna, bamRna)
        sub_checked_read_presence <- lapply(sub_checked_read_presence_list,
                                            nth, n=1) %>%
          bind_rows() %>%
          filter(nreads>0)
        bam_raw <- lapply(sub_checked_read_presence_list,
                          nth, n=2)
        ## now check with the remaining ones if there is a connection
        ## first check if the main variants are still there
        still_present_muts <- c(sub_checked_read_presence$mut_id1,
                                sub_checked_read_presence$mut_id2) %>%
          unique()
        if(main_comb[["mut_id1"]] %in% still_present_muts&
           main_comb[["mut_id2"]] %in% still_present_muts){
          
          ## make list to use with A* algorithm
          path <- find_path(still_present_muts, sub_checked_read_presence,
                            main_comb)
          if(!is.null(path)){
            catt(printLog, 5, c("path found:", 
                paste(unlist(path), collapse=" - ")))
            ## find relevant combinations to phase
            if(printLog==TRUE){
              apply(
                vars_in_between[
                  which(vars_in_between$mut_id %in% unlist(path)),],
                1, function(r){
                  catt(printLog, 5, c(r[["mut_id"]],
                                      r[["chr"]],
                                      r[["pos"]]))
                })
            }
            sub_fin_comb <- phase_along_path(path, sub_checked_read_presence,
                                             bam_raw, purity, printLog)
            non_null <- sub_fin_comb %>%
              filter(status!="null")
            if(nrow(non_null)!=0&
               (main_comb[["mut_id1"]] %in% c(non_null$f, non_null$l)&
                main_comb[["mut_id2"]] %in% c(non_null$f, non_null$l))){
              catt(printLog, 5, 
                   "all combinations of path have assigned status")
              final <- recombine_phasing_results(non_null, main_comb)
              ext_snp_phasing <- sub_fin_comb %>%
                mutate(main_comb=main_comb[['comb_id']])
              if(nrow(final)!=0){
                return_null_result <- FALSE
                RESULT <- finalize_snp_phasing(final, main_comb)
                full_read_info <- NULL    
                ext_snp_info <- paste(
                  "SNP phasing successfull:", RESULT$status, "; wt cp:",
                  RESULT$wt_cp)
              } else {
                ext_snp_info <- 
                  "recombination impossible: more than one diff status in path"
              }
            } else {
              ext_snp_info <- "recombination imposiible: null status in path"
            }
          } else {
            ext_snp_info <- "no connecting path between main muts"
          }
        } else {
          ext_snp_info <- "no overlap with main muts"
        }
      } else {
        ext_snp_info <- "no SNPs between"
      }
    } else { 
      ext_snp_info <- "muts too far away for SNP phasing"
    }
    catt(printLog, 3, ext_snp_info)
    if(return_null_result==TRUE){
      RESULT <- return_null_result(main_comb, bamRna, ext_snp_info)
      full_read_info <- NULL
    }
    return(list(assigned=RESULT, detailed=full_read_info, 
                snp_phasing=ext_snp_phasing))
  }) 
  return(final_combinations) 
}
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_replace_all str_count str_detect
#' @importFrom dplyr as_tibble group_by relocate bind_rows left_join select tally mutate filter nth relocate
predict_zygosity_genewise <- function(GENE, df_all_mutations, bamDna, 
                                      bamRna, showReadDetail,
                                      printLog, purity, vcf, distCutOff){
  if(printLog==TRUE){
    message(GENE)
  }
  gene <- pre_info <-  chr <-  pos <- wt_cp <-  mut_id <- status <- . <- 
    score <- comb <- dist <- tcn <- info <- n <- NULL 
  start_gene_eval <- Sys.time()
  pre_df_gene <- df_all_mutations %>% filter(gene==GENE)    
  df_gene_raw <- pre_df_gene%>%
    # all germline variants which are lost in the tumor are excluded
    filter(!str_detect(pre_info, 'variant lost in tumor'))
  ## check if variants at the same position are present
  n_vars_per_pos <- df_gene_raw %>%
    group_by(chr, pos) %>%
    tally() %>%
    filter(n>1)
  if(nrow(n_vars_per_pos)!=0){
    var_to_remove <- apply(n_vars_per_pos, 1, function(VAR){
      selected <- df_gene_raw %>% 
        filter(chr==VAR[["chr"]],
               pos==VAR[["pos"]]) %>%
        filter(wt_cp==min(.$wt_cp)) %>%
        ## if both vars have the same impac we just select the first
        .[1,]
      warning(
          "more than one variant detected at position:", 
          VAR[["chr"]], ":", VAR[["pos"]],
          "\n  maybe an indel with shifted annotation of position?",
          "\nSelecting variant of highest relevance by lowest left wt-copies:",
          selected$mut_id)
      to_remove <- df_gene_raw %>% 
        filter(chr==VAR[["chr"]],
               pos==VAR[["pos"]]) %>%
        filter(mut_id!=selected$mut_id) %>%
        select(chr, pos, mut_id) %>%
        mutate(remove=TRUE)
      return(to_remove)
    }) %>%
      bind_rows()
    df_gene <- df_gene_raw %>%
      left_join(var_to_remove, 
                by=c("chr"="chr", "pos"="pos", "mut_id"="mut_id")) %>%
      filter(is.na(remove)) %>%
      select(-remove)
  } else {
    df_gene <- df_gene_raw
  }
  if(!nrow(df_gene)==0){
    log_obj <- list(gene=GENE, n_mut=nrow(df_gene))
    all_comb <- NULL
    all_comb_to_export <- NULL
    read_details <- NULL
    snp_phasing <- NULL
## at first check if any of the mutations already affects all alleles           
    if(str_detect(paste(df_gene$pre_info, collapse = ' '), 
                  'all copies affected')){
      concern_info <- df_gene %>% 
        filter(str_detect(pre_info, "all copies affected")) %>% .[1,]
      class <- get_classification(concern_info)
      origin <- ifelse(str_detect(concern_info$origin, "germ"),
                       "germ-",
                       "som-")
      log_obj[['score']] <- 2
      log_obj[['info']] <- concern_info$pre_info                 
      log_obj[['tcn']] <- concern_info$tcn
      log_obj[['aff_cp']] <- as.numeric(concern_info$aff_cp) 
      log_obj[['class']] <- paste0(origin, class)
      catt(printLog, 1, c("all copies affected by variant:", 
            concern_info$mut_id))
      catt(printLog, 1, concern_info$pre_info)
    } else if(nrow(df_gene)==1){
      log_obj[['score']] <- 1
      log_obj[['info']] <- df_gene$pre_info
      log_obj[['tcn']] <- df_gene$tcn
      log_obj[['aff_cp']] <- as.numeric(df_gene$aff_cp)
      log_obj[['class']] <- paste0(get_classification(df_gene))
      catt(printLog, 1, "1 variant detected:")
      catt(printLog, 1, df_gene$pre_info)
    } else {  
      catt(printLog, 1, c(nrow(df_gene) ,
            "variants detected: Initializing haplotype phasing"))
      all_comb_raw <- phase(df_gene, bamDna, bamRna, showReadDetail,
                            purity, vcf, distCutOff, printLog)
      all_comb <- lapply(all_comb_raw, nth, n=1) %>% 
        Reduce(function(x,y)rbind(x,y),.) %>%
        mutate(gene=GENE)
      all_comb_to_export <- all_comb %>% 
        select(-score) %>%
        relocate(gene, comb, dist, tcn, 
                        status, left_wt_cp=wt_cp)
      snp_phasing <- lapply(all_comb_raw, nth, n=3) %>%   
        bind_rows() %>%
        mutate(gene=GENE)
      if(showReadDetail==TRUE){
        read_details <- lapply(all_comb_raw, nth, n=2) %>%   
          bind_rows() %>%
          mutate(gene=GENE)        
      }
      if(nrow(all_comb)==1){
        log_obj[['score']] <- all_comb$score
        log_obj[['info']] <- all_comb$info
        log_obj[['tcn']] <- all_comb$tcn
        log_obj[['aff_cp']] <- as.numeric(all_comb$tcn-all_comb$wt_cp)
        log_obj[['class']] <- paste("comb:", all_comb$comb)
      } else if(max(all_comb$score)==2){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 2
        log_obj[['info']] <- 
          "at least one of the combinations of muts affects all copies"
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(str_count(paste(all_comb$status, collapse=' '),
                          'diff')==nrow(all_comb)&
                nrow(all_comb)>=mean(as.numeric(df_gene$tcn))){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- paste(
          'unclear: all affected is very likely because all',nrow(all_comb),
          'mutations are at different copies and tcn=', 
          mean(as.numeric(df_gene$tcn)),
          'but as said in the meeting we are not sure if',
          'tumor is subclonal and therefore wt-copies left')
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(str_detect(paste(all_comb$status, collapse=' '),
                           'unclear: muts at different copies and tcn=2')){
        most_important_comb <- all_comb %>% 
          filter(str_detect(status, 
                            'unclear: muts at different copies and tcn=2')) %>%
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- paste(
          "unclear: all affected is very likely because two mutations were",
          "detected at different reads and the tcn is 2 but", 
          "we are not sure if tumor is subclonal and therefore only ", 
          "wt-copies left")
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else if(nrow(all_comb)>3){
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- 
          paste('unclear: very high number of mutations...', 
                'probably wt copies left but please check variant table')
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", most_important_comb$comb)
      } else {
        most_important_comb <- all_comb %>% 
          filter(wt_cp==min(all_comb$wt_cp)) %>%
          filter(rownames(.)==1)
        log_obj[['score']] <- 1
        log_obj[['info']] <- 'wt-copies left' 
        log_obj[['tcn']] <- most_important_comb$tcn
        log_obj[['aff_cp']] <- 
          most_important_comb$tcn-most_important_comb$wt_cp %>% 
          as.numeric()
        log_obj[['class']] <- paste("comb:", 
                                    most_important_comb$comb)            
      }
      log_obj[['info']] <- paste("phasing:", log_obj[['info']])
    } 
    df_reduced <- log_obj %>% as_tibble() %>%
      mutate(status=ifelse(score==2,
                           "all_copies_affected",
                           "wt_copies_left"),
             info=str_replace_all(info, 
                                  " -> all copies affected| -> wt copies left",
                                  "")) 
    zygosity_gene <- list(pre_df_gene, 
                          list(df_reduced, 
                               all_comb_to_export), 
                          read_details, snp_phasing)
  } else {
    zygosity_gene <- list(pre_df_gene, NULL)
  } 
  if(printLog==TRUE){
    end_gene_eval <- Sys.time()
    message("~", round(
                    as.numeric(
                      difftime(
                        end_gene_eval, 
                        start_gene_eval,
                        units="secs")
                      ),3), "s") 
  }
  return(zygosity_gene)
}
#' @keywords internal
#' description follows 
#' @importFrom dplyr case_when
get_classification <- function(data){
  case_when(
    data$class=="homdel" ~ "homdel",
    nchar(data$alt)==nchar(data$ref) ~ "snv",
    nchar(data$alt)>nchar(data$ref) ~ "ins",
    nchar(data$alt)<nchar(data$ref) ~ "del",
  ) %>%
    return()
}
#' @keywords internal
#' description follows
#' @importFrom dplyr case_when
define_class <- function(ref, alt){
  case_when(
    nchar(alt)==nchar(ref) ~ "snv",
    nchar(alt)>nchar(ref) ~ "ins",
    nchar(alt)<nchar(ref) ~ "del",
  ) %>%
    return()
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr as_tibble filter select mutate mutate_all
combine_uncovered_input_variants <- function(somSmallVars, germSmallVars,
                                             som_covered, germ_covered,
                                             templateGenes){
  mid <- seqnames <- start <- ref <- alt <- gene <- uncovered_som <- 
    uncovered_germ <- NULL
  if(!is.null(somSmallVars)){
    uncovered_som <-  as_tibble(somSmallVars) %>%
      filter(!mid %in% som_covered) %>%
      select(chr=seqnames, pos=start,ref, alt, gene) %>%
      mutate(origin="somatic") %>%
      mutate_all(.funs = as.character)
  } 
  if(!is.null(germSmallVars)){
    uncovered_germ <- as_tibble(germSmallVars) %>%
      filter(!mid %in% germ_covered) %>%
      select(chr=seqnames, pos=start,ref, alt, gene) %>%
      mutate(origin="germline") %>%
      mutate_all(.funs = as.character)
  }
  
  if(nrow(bind_rows(uncovered_som, uncovered_germ))==0){
    combined_uncovered <- NULL
  } else {
    combined_uncovered <- bind_rows(uncovered_som, uncovered_germ) %>%
      filter(gene %in% templateGenes)
  }
  return(combined_uncovered)
}
#' @keywords internal
check_opt_assgap <- function(assumeSomCnaGaps, ploidy){
  if(assumeSomCnaGaps==TRUE&is.null(ploidy)){
    warning("somatic CNA gaps can only be assumed if input ploidy is",
            "provided. Provide ploidy=2 to assume diploid case")
    return(FALSE)
  } else {
    return(assumeSomCnaGaps)
  }
}
#' @keywords internal
check_opt_incdel <- function(includeIncompleteDel, ploidy){
  if(includeIncompleteDel==TRUE&is.null(ploidy)){
    warning("Large scale deletions cannot be included without ploidy",
                "Please provide input ploidy. Provide ploidy=2",
                "to assume diploid case")
    return(FALSE)
  } else {
    return(includeIncompleteDel)
  }
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr as_tibble group_by mutate ungroup filter
#' @importFrom purrr compact
combine_main_variant_tables <- function(df_germ, df_som, df_homdels,
                                        templateGenes, purity){
  gene <- . <- NULL
  df_all_mutations_unfiltered <- list(df_germ, df_som, df_homdels) %>% 
    compact() %>%
    Reduce(function(x,y)rbind(x,y),.) %>% 
    as_tibble() %>% 
    group_by(gene) %>% 
    mutate(mut_id=paste0("m", seq_len(length(gene))))%>% 
    ungroup() %>%
    mutate(purity=purity) 
  df_all_mutations <- df_all_mutations_unfiltered %>%
    filter(gene %in% templateGenes)
  if(nrow(df_all_mutations_unfiltered)!=nrow(df_all_mutations)){
    warning(
      abs(nrow(df_all_mutations_unfiltered)-nrow(df_all_mutations)),
      "variants in genes not found in reference genome annotation",
      "(geneModel). They are removed from the analysis.",
      "Please provide a annotation that contains every gene from the",
      "variant inputs"
    )
  }
  return(df_all_mutations)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr bind_rows group_by mutate relocate select ungroup relocate
bind_incdel_to_pre_eval <- function(df_incompletedels, df_all_mutations){
  purity <- gene <- mut_id <- NULL
  if(is.null(df_all_mutations)){
    full_df_all_mutations <- NULL
  } else {
    if(is.null(df_incompletedels)){
      full_df_all_mutations <- df_all_mutations %>%
        select(-purity) %>%
        relocate(gene, mut_id, 2:13)
    } else {
      full_df_all_mutations <- bind_rows(
        df_all_mutations,
        df_incompletedels %>% mutate(mut_id=NA)
      ) %>% group_by(gene) %>% 
        mutate(mut_id=paste0("m", seq_len(length(gene))))%>% 
        ungroup()%>%
        select(-purity) %>%
        relocate(gene, mut_id, 2:13)       
    }    
  }

  return(full_df_all_mutations)
}
#' @keywords internal
#' @importFrom stringr %>%
#' @importFrom dplyr filter bind_rows mutate select
bind_incdel_to_final_eval <- function(df_incompletedels, final_output){
  gene <- 0
  if(is.null(final_output)){
    full_output <- NULL
  } else {
    if(is.null(df_incompletedels)){
      full_output <- final_output
    } else {
      full_output <- bind_rows(
        final_output,
        df_incompletedels %>%
          filter(!gene %in% final_output$gene) %>%
          select(gene) %>%
          mutate(status="wt_copies_left",
                 info="somatic-incompletedel")
      )       
    }    
  }
  return(full_output)
}
#' predicts zygosity of a set of genes of a sample
#' @param purity purity of the sample (numeric value between 0 and 1 indicating 
#' the fraction of relevant sample with control/unrelevant tissue)
#' @param ploidy ploidy of the sample (numeric value)
#' @param sex sex of the sample (character: "male", "female", "m", "f")
#' @param somCna GRanges object containing all genomic regions with annotated 
#' total copynumber and cna_type as metadata columns. The total-copynumber 
#' column should be named "tcn" but also some other commonly used names. 
#' It should contain numeric values or characters that can be converted to 
#' numeric values. The cna_type column must contain the information about 
#' loss of heterozygosity (LOH). Therefore the term "LOH" must be explicitely 
#' mentioned in the column. If a genomic region is not present in the object, 
#' it will be taken as heterozygous with neutral TCN of 2. 
#' @param somSmallVars GRanges object containing all somatic small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), 
#' alternative base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF). 
#' If the object is not provided the tool assumes there are no somatic small 
#' variants.
#' @param germSmallVars GRanges object containing all germline small 
#' variants (SNV and INDEL).
#' Required metadata columns are reference base (ref/REF), alternative 
#' base (alt/ALT),
#' annotation of the gene name (gene/GENE) and the allele-frequency (af/AF)
#' If the object is not provided the tool assumes there are no germline small 
#' variants.
#' @param geneModel GRanges object containing the gene-annoattion of 
#' the used reference genome with metadata column of the gene name (gene)
#' @param bamDna path to bam-file
#' @param bamRna optional; path to rna file (bam format)
#' @param includeHomoDel default = TRUE; if FALSE homozygous deleteions are 
#' excluded
#' @param includeIncompleteDel default = TRUE; if FALSE heterzygous deleteions 
#' are excluded
#' @param showReadDetail default = FALSE; if TRUE a table is added to the 
#' output,
#' @param printLog default = FALSE; if TRUE the gene which is evaluated is 
#' printed in console, 
#' containing the query-name of each read which was used to perform 
#' haplotype-phasing and the info into which class it was assigned.
#' @param assumeSomCnaGaps (logical, default=FALSE) Only required if the somCna
#' object lacks copy number information for genomic segments on which small 
#' variants are detected. By default, variants in such regions will be excluded 
#' from the analysis as required information about the copy number is missing. 
#' These variants will be attached to the final output list in a separate 
#' tibble. To include them, this flag must be set TRUE and the ground ploidy 
#' must be given as an input. This ground ploidy will then be taken as tcn in 
#' the missing regions. If no ploidy is given the tool will assume the ground 
#' ploidy of 2 when this flag is TRUE.
#' @param byTcn logical, default=TRUE; optional if includeHomoDel or 
#' includeIncompleteDelS is TRUE. If FALSE the tool will not use tcn as a 
#' criterion to assign large deletions. It will use the cna_type column and 
#' check for indicating strings like HOMDEL/HomoDel/DEL. Some commonly used 
#' strings are covered. It is recommended to leave this flag TRUE
#' @param colnameTcn character indicating the name of the metadata containing 
#' the tcn information in the somCna object. If not provided the tool tries to 
#' detect the column according to default names
#' @param colnameCnaType character indicating the name of the metadata 
#' containing cna type information in the somCna object. 
#' If not provided the tool tries to detect the column according to default 
#' names
#' @param vcf character; path to variant call file (.vcf.gz format). 
#' Will be used (if provided)
#' for extended SNP phasing if variants on the same gene are too far away from
#' each other for direct haplotype phasing
#' @param distCutOff numeric, default=5000; if input vcf is provided and SNP
#' phasing is performed, this will limt the distance at which the SNP phasing
#' should not be tried anymore. As the probability of finding overlapping reads
#' at such a long distance is very low and the runtime will increase
#' exponentially.
#' @return A list of dataframes. Those are the evaluation per variant, 
#' the evaluation per gene and, if performed, the info about the 
#' haplotype-phasing.
#' @examples
#' cnvs  = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr = "chr17",
#'     start = c(170060, 34520990),
#'     end = c(34520990, 83198614),
#'     tcn = c(2, 1),
#'     cna_type = c("neutral", "LOH")
#'   )
#' )
#' somatic_vars = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr="chr17",
#'     start = 7675088,
#'     end = 7675088,
#'     ref = "C",
#'     alt = "T",
#'     af = 0.65,
#'     gene = "TP53" 
#'   )
#' )
#' germline_vars = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr="chr17",
#'     start = 41771694,
#'     end = 41771694,
#'     ref = "GTGT",
#'     alt = "G",
#'     af = 0.95,
#'     gene = "JUP" 
#'   )
#' )
#' reference = GenomicRanges::GRanges(
#'   dplyr::tibble(
#'     chr = "chr17",
#'     start = c(7661778, 41754603),
#'     end = c(7687538, 41786931),
#'     gene = c("TP53", "JUP")
#'   )
#' )
#' sex = "female"
#' purity = 0.9
#' bamfile <- system.file("extdata", "ZP_example.bam", 
#'   package = "ZygosityPredictor")
#' predict_zygosity(purity = purity, sex = sex, 
#'   somCna = cnvs,
#'   somSmallVars = somatic_vars,
#'   germSmallVars = germline_vars,
#'   geneModel = reference,
#'   bamDna = bamfile
#' )
#' @importFrom stringr %>%
#' @importFrom IRanges subsetByOverlaps
#' @importFrom purrr compact
#' @importFrom dplyr bind_rows nth select tibble
#' @export
predict_zygosity <- function(purity, 
                             sex,
                             somCna, 
                             geneModel,
                             bamDna,
                             somSmallVars=NULL, 
                             germSmallVars=NULL,
                             bamRna=NULL,
                             ploidy=NULL, 
                             colnameTcn=NULL,
                             colnameCnaType=NULL,
                             includeHomoDel=TRUE,
                             includeIncompleteDel=TRUE,
                             showReadDetail=FALSE,
                             printLog=FALSE,
                             assumeSomCnaGaps=FALSE,
                             byTcn=TRUE,
                             vcf=NULL,
                             distCutOff=5000
                             ){
  status <- info <- wt_cp <- . <- df_homdels <- df_all_mutations <- gene <-
    final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- 
    combined_snp_phasing <- NULL
  ## check input for valid format
  somSmallVars <- check_gr_small_vars(somSmallVars, "somatic")
  germSmallVars <- check_gr_small_vars(germSmallVars, "germline")
  purity <- check_purity(purity)
  ploidy <- check_ploidy(ploidy)
  sex <- check_sex(sex)
  bamDna <- check_bam(bamDna)
  geneModel <- check_gr_gene_model(geneModel)
  assumeSomCnaGaps <- check_opt_assgap(assumeSomCnaGaps, ploidy)
  somCna <- check_somCna(somCna, geneModel, sex, ploidy, 
                         assumeSomCnaGaps, colnameTcn, 
                         colnameCnaType) 
  includeIncompleteDel <- check_opt_incdel(includeIncompleteDel, ploidy)
  includeHomoDel <- check_opt_incdel(includeHomoDel, ploidy)
  bamRna <- check_rna(bamRna)
  vcf <- check_vcf(vcf)
  templateGenes <- geneModel$gene
  ## get variants not covered by CNV input
  if(!is.null(somSmallVars)){
    gr_som_cov <- subsetByOverlaps(somSmallVars, somCna) 
    som_covered <- gr_som_cov$mid
  }
  if(!is.null(germSmallVars)){
    gr_germ_cov <- subsetByOverlaps(germSmallVars, somCna) 
    germ_covered <- gr_germ_cov$mid
  }
  combined_uncovered <- combine_uncovered_input_variants(somSmallVars, 
                                                         germSmallVars,
                                                         som_covered,
                                                         germ_covered,
                                                         templateGenes)
  ## preapre input data for prediction (pre-evaluation)
  df_germ <- prepare_germline_variants(
    gr_germ_cov, 
    somCna, purity, sex)
  df_som <- prepare_somatic_variant_table(
    gr_som_cov, 
    templateGenes, 
    somCna, purity, sex)
  df_homdels <- extract_all_dels_of_sample(somCna, geneModel, 
                                             "homdel", byTcn, sex, ploidy,
                                             includeHomoDel)
  ## predict zygosity for each gene
  if(!is.null(df_germ)|!is.null(df_som)|!is.null(df_homdels)){
    df_all_mutations <- combine_main_variant_tables(df_germ, df_som, 
                                                    df_homdels,
                                                    templateGenes, purity)
    result_list <- lapply(
        unique(df_all_mutations$gene), 
        predict_zygosity_genewise, 
        df_all_mutations, 
        bamDna,
        bamRna,
        showReadDetail,
        printLog,
        purity,
        vcf,
        distCutOff)
    pre_scoring <- lapply(result_list, nth, n=1) %>% 
      bind_rows()
    final_result_list <- lapply(result_list, nth, n=2) %>% compact()
    if(length(final_result_list)!=0){
      final_phasing_info <- lapply(final_result_list, nth, n=2) %>% 
        compact() %>% 
        bind_rows() 
      final_output <- lapply(final_result_list, nth, n=1) %>% 
        bind_rows() %>%
        select(gene, status, info)
      combined_read_details <- lapply(result_list, nth, n=3) %>% compact() %>%
        bind_rows()
      combined_snp_phasing <- lapply(result_list, nth, n=4) %>% compact() %>%
        bind_rows()
    } 
  } 
  ## include incomplete dels to output
  df_incompletedels <- extract_all_dels_of_sample(somCna, geneModel, 
                                              "incompletedel", TRUE, sex, 
                                              ploidy, includeIncompleteDel) 
  evaluation_per_variant <- bind_incdel_to_pre_eval(df_incompletedels,
                                                   df_all_mutations)
  evaluation_per_gene <- bind_incdel_to_final_eval(df_incompletedels,
                                                   final_output)
  ## export output
  result_list <- list(
    eval_per_variant=evaluation_per_variant,
    eval_per_gene=evaluation_per_gene,
    phasing_info=final_phasing_info,
    readpair_info=combined_read_details,
    uncovered_input=combined_uncovered,
    ext_snp_phasing=combined_snp_phasing
  ) %>%
    compact()
  return(result_list)
} 
## the following functions were taken from 
## https://github.com/machow/astar-r
#' @keywords internal
make_search_node <- function(data, gscore, fscore) {
  env <- new.env()
  env$data <- data
  env$gscore <- gscore
  env$fscore <- fscore
  env$closed <- FALSE
  env$out_openset <- TRUE
  env$came_from <- NULL
  env
}
#' @keywords internal
reconstruct_path <- function(goal) {
  path <- list(goal$data)
  crnt <- goal
  while (!is.null(crnt$came_from)) {
    crnt <- crnt$came_from
    path <- c(list(crnt$data), path)
  }
  path
}
#' @keywords internal
#' @importFrom datastructures binomial_heap insert peek pop
a_star_pathfinder <- function(nodes, start, goal,
                              hash_func = identity, search_node_env = NULL) {
  if(start==goal)
    return(list(start))
  search_nodes <- if (!is.null(search_node_env)) search_node_env else list()
  start_node <- make_search_node(start, gscore = 0,
                                 fscore = 1)
  start_hash <- hash_func(start)
  search_nodes[[start_hash]] <- start_node
  open_set <- binomial_heap("numeric")
  insert(open_set, start_node$fscore, start_hash)
  while (!is.null(peek(open_set))) {
    crnt <- search_nodes[[pop(open_set)[[1]]]]
    #if (identical(crnt$data, goal))
    if(crnt$data==goal)
      return(reconstruct_path(crnt))
    crnt$out_openset <- TRUE
    crnt$closed <- TRUE
    for (neighbor in names(nodes[[crnt$data]])) {
      indx <- hash_func(neighbor)
      neigh_node <- search_nodes[[indx]]
      if (is.null(neigh_node)) {
        neigh_node <- search_nodes[[indx]] <-
          make_search_node(neighbor, Inf, Inf)
      }
      if (neigh_node$closed) next
      tentative_gscore <- crnt$gscore + nodes[[crnt$data]][neigh_node$data]
      if (tentative_gscore >= neigh_node$gscore) next
      neigh_node$came_from <- crnt
      neigh_node$gscore <- tentative_gscore
      neigh_node$fscore <-
        tentative_gscore + 1
      if (neigh_node$out_openset) {
        neigh_node$out_openset <- FALSE
        insert(open_set, neigh_node$fscore, indx)
      }
    }
  }
}