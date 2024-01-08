#' @keywords internal
#' description follows
#' @importFrom stringr str_match
#' @importFrom dplyr case_when
#' @importFrom GenomicRanges elementMetadata elementMetadata<-
assign_correct_colnames <- function(obj, type){
  func_start()
  #vm(as.character(sys.call()[1]), 1)
  . <- NULL
  if(type=="scna"){
    
    ## define correct tcn column
    col_tcn <- str_match(
      nm_md(obj), paste(allowed_inputs("colnames_tcn", regex=TRUE), collapse="|")) %>%
      .[which(!is.na(.))]
    elementMetadata(obj)[,"tcn"] <- 
      elementMetadata(obj)[,col_tcn] 
    ## define correct cna_type column
    if("LOH" %in% nm_md(obj)){
      obj$cna_type <- case_when(obj$LOH==TRUE ~ "LOH",
                                TRUE ~ "HZ")
    } else {
      col_cna_type <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_cna_type", regex=TRUE), 
                          collapse="|")) %>%
        .[which(!is.na(.))]   
      elementMetadata(obj)[,"cna_type"] <- 
        elementMetadata(obj)[,col_cna_type]
      obj$LOH <- case_when(str_detect(obj$cna_type, "LOH") ~ TRUE,
                           TRUE ~ FALSE)
    }
    ## define allelic imbalance column (optional input)
    ## check if column provided in input
    col_all_imb <- str_match(
      nm_md(obj), paste(allowed_inputs("colnames_all_imb", regex=TRUE), 
                        collapse="|")) %>%
      .[which(!is.na(.))] 
    if(length(col_all_imb)>0){
    ## if colum already there, rename it
      elementMetadata(obj)[,"all_imb"] <- 
        elementMetadata(obj)[,col_all_imb]
    } else {
    ## if not yet there, calculate it for segments where it is possible
      obj$all_imb <- case_when(
        obj$LOH==FALSE&round(as.numeric(obj$tcn))%%2!=0&as.numeric(obj$tcn)>=2.5 ~ TRUE,
        TRUE ~ FALSE
      )
      #obj$gt_cna <- "1:2"
    }
    ## if allelic imbalance was detected somewhere, also the genotype needs to be annotated 
    ## except if tcn ==3
    col_genotype <- str_match(
      nm_md(obj), paste(allowed_inputs("colnames_genotype", regex=TRUE), 
                        collapse="|")) %>%
      .[which(!is.na(.))]
    if(length(col_genotype)>0){
      ## if colum already there, rename it
      ## and check if it has only entries that are actually imbalanced... remove
      ## something like 2:2 or 3:3
      input_gt_col <- as.character(elementMetadata(obj)[,col_genotype]) %>%
        lapply(.,function(GT){
          uniqued_gts <- str_split(GT, ":") %>% unlist() %>% unique()
          if(length(uniqued_gts)==1){
            warning("column gt_cna contains annotations of balanced segments. They are removed from imbalance phasing if enabled")
            return(NA)
          } else {
            return(GT)
          }
        }) %>%
        unlist()
      elementMetadata(obj)[,"gt_cna"] <- input_gt_col
    } else {
      ## if not yet there, calculate it for segments where it is possible
      obj$gt_cna <- case_when(
        round(obj$tcn)==3 ~ "1:2",
        TRUE ~ NA
      )
    }
  } else {
    col_gene <- str_match(
      nm_md(obj), 
      paste(allowed_inputs("colnames_gene", regex=TRUE), collapse="|")) %>%
      .[which(!is.na(.))]
    elementMetadata(obj)[,"gene"] <- 
      elementMetadata(obj)[,col_gene]  
    if(type=="small_vars"){
      col_af <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_af", regex=TRUE), collapse="|")) %>%
        .[which(!is.na(.))]
      col_ref <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_ref", regex=TRUE), collapse="|")) %>%
        .[which(!is.na(.))]
      col_alt <- str_match(
        nm_md(obj), paste(allowed_inputs("colnames_alt", regex=TRUE), collapse="|")) %>%
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
  }
  func_end()
  return(obj)
}
#' @keywords internal
#' no dependencies
formula_checks <- function(chr, af, tcn, purity, sex, c_normal, af_normal=0.5){
  purity <- check_purity(purity)
  af <- check_af(af)
  af_normal <- check_af(af_normal)
  tcn <- check_tcn(tcn)
  if(is.null(c_normal)){
    sex <- check_sex(sex)
    chr <- check_chr(chr)
    if((sex=="male"&(chr=="X"|chr=="chrX"))|(chr=="Y"|chr=="chrY")){
      c_normal <- 1
      af_normal <- 1
    } else {
      c_normal <- 2
      af_normal <- 0.5
    }
  } else {
    c_normal <- check_ploidy(c_normal)
  }
  return(list(af=af, tcn=tcn, purity=purity, c_normal=c_normal, 
              af_normal=af_normal))
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
#' returns names of metadata columns from GRanges object
#' @importFrom GenomicRanges elementMetadata
nm_md <- function(obj){
  return(names(elementMetadata(obj)))
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
#' @keywords internal
#' description follows
#' @importFrom stringr %>% str_detect
#' @importFrom GenomicRanges elementMetadata elementMetadata<- seqnames
#' @importFrom methods is
check_somCna <- function(somCna, geneModel, sex, ploidy,
                         assumeSomCnaGaps, colnameTcn, 
                         colnameCnaType){
  func_start()
  . <- NULL
  somCna <- general_gr_checks(somCna, "scna", "somCna")
  somCna$tcn_assumed <- FALSE
  ## check if there are segments without cna_type annotation
  if(sum(is.na(elementMetadata(somCna)[,"cna_type"]))>0){
    warning("cna_type column of input somCna contains", 
            sum(is.na(elementMetadata(somCna)[,"cna_type"])),
            "NA values;\n  they will be taken as hetero-zygous",
            "/hemizygous for gonosomes in male samples\n")
    elementMetadata(somCna)[,"cna_type"][
      which(
        is.na(elementMetadata(somCna)[,"cna_type"]))] <- NA
  }
  ## for male samples the gonosomes are hemizygous, we treat them like autosomal
  ## region with LOH, therefore we annotate LOH to the cna_type column
  if(sex=="male"&str_detect(
    paste(as.character(seqnames(somCna)), collapse=" "), "X|Y")){
    ## if true, the sample is male and has Gonosomal regions
    somCna[which(
      as.character(seqnames(somCna)) %in% c("X", "Y"))]$cna_type <-
      paste0(
        somCna[which(
          as.character(seqnames(somCna)) %in% c("X", "Y"))]$cna_type,
        ";LOH")
  }
  ## if CNA gaps should be assumed to ploidy
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
  ## if cna gaps should not be assumed, missing tcn segments will be excluded
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
  new_somCna$seg_id <- seq(1, length(new_somCna))
  func_end()
  return(new_somCna)
}
#' @keywords internal
#' description follows
check_name_presence <- function(obj, type){
  func_start()
  if(type=="scna"){
    if(
      any(allowed_inputs("colnames_tcn") %in% nm_md(obj))&
      (any(allowed_inputs("colnames_cna_type") %in% nm_md(obj))|
       "LOH" %in% nm_md(obj))
    ){ 
      res <- TRUE
    } else {
      res <- FALSE
    } 
  } else {
    if(
      (any(allowed_inputs("colnames_gene") %in% nm_md(obj))&
        (type=="gene_model"|
         (
           any(allowed_inputs("colnames_af") %in% nm_md(obj))&
           any(allowed_inputs("colnames_ref") %in% nm_md(obj))&
           any(allowed_inputs("colnames_alt") %in% nm_md(obj))
         )
        )
      )
    ){
      res <- TRUE
    } else {
      res <- FALSE
    }    
  }
  func_end()
  return(res)
}
check_logDir <- function(logDir){
  if(!is.null(logDir)){
    if(file.exists(logDir)){
      return(logDir)
    } else {
      warning("logDir does not exist, no log will be stored")
      return(NULL)
    }
  } else {
    return(NULL)
  }
}
general_gr_checks <- function(obj, type, lab){
  func_start()
  if(!is(obj, "GRanges")){
    stop("input ", lab, " must be a GRanges object;",
         "given input appears to be:", 
         class(obj))
  } else if(type=="haploBlocks"){
  ## haploblocks do not need any metadata columns  
    res <- obj
  } else if(!check_name_presence(obj, type)){
    if(type=="scna"){
      stop(
        "input somCna requires the following metadata columns: ",
        "\'tcn\' and \'cna_type\'")
    } else if(type=="small_vars"){
      stop("input ", lab, " requires the following metadata columns:",
           "\'gene\'/\'GENE\', \'ref\'/\'REF\',",
           " \'alt\'/\'ALT\' and \'af\'/\'AF\'")      
    } else {
      stop(
        "input geneModel requires the following metadata columns: ",
        "\'gene\'/\'GENE\'")
    }
  } else {
    res <- assign_correct_colnames(obj, type)
  }
  func_end()
  return(res)
}
#' @keywords internal
#' description follows
check_gr_gene_model <- function(geneModel, is_pre_eval){
  func_start()
  . <- NULL
  if(is_pre_eval==TRUE){
    if(is.null(geneModel)){
      res <- NULL
    } else {
      res <- general_gr_checks(geneModel, "gene_model", "geneModel")
    }
  } else if(is.null(geneModel)){
    stop("input geneModel must not be NULL")
  } else {
    res <- general_gr_checks(geneModel, "gene_model", "geneModel")
  }
  func_end()
  return(res)
}
#' @keywords internal
#' description follows
check_gr_small_vars <- function(obj, origin){
  func_start()
  . <- NULL
  lab <- ifelse(origin=="somatic",
                "somSmallVars",
                "germSmallVars")
  if(is.null(obj)){
    warning("Input ", lab, " empty/does not contain variants. ",
            "Assuming there are no ", origin, " small variants")
    res <- NULL
  } else if(length(obj)==0){
    warning("Input ", lab, " empty/does not contain variants. ",
            "Assuming there are no ", origin, " small variants")
    res <- NULL
  } else {
    res <- general_gr_checks(obj, "small_vars", lab)
  }
  func_end()
  return(res)
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
  func_start()
  if(is.null(ploidy)){
    res <- NULL
  } else {
    if(is.na(as.numeric(ploidy))){
      stop("input ploidy/c_normal must be numeric or a character that can",
           "be converted to numeric;\n  ", ploidy, 
           "can not be converted to numeric")
    } else {
      res <- as.numeric(ploidy)
    }
  }
  func_end()
  return(res)
}
#' @keywords internal
#' description follows
#' @importFrom stringr str_detect
check_sex <- function(sex){
  . <- NULL
  allowed_sex <- allowed_inputs("sex_names")
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
#' @importFrom Rsamtools TabixFile
check_vcf <- function(vcf){
  if(is.null(vcf)){
    return(NULL)
  } else if(!sum(unlist(lapply(vcf, file.exists)))==length(vcf)){
    stop("one of input vcfs does not exist")
  } else {
    processed_vcf <- lapply(vcf, function(VCF){
      if(file.exists(paste0(VCF, ".tbi"))){
        return(Rsamtools::TabixFile(VCF))
      } else {
        return(VCF)
      }
    })
    return(processed_vcf)
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
check_haploblocks <- function(haploBlocks){
  if(is.null(haploBlocks)){
    return(NULL)
  } else {
    haploBlocks <- general_gr_checks(haploBlocks, "haploBlocks", "haploBlocks")
    haploBlocks$hap_id <- seq(1, length(haploBlocks))
    return(haploBlocks)    
  }
}
#' @keywords internal
allowed_inputs <- function(which_one, regex=FALSE){
  . <- NULL
  type_list <- list(
    cna_homdel_annotation = paste(c("HOMDEL", "HomoDel", "HomoDEL", "HOMODEL", 
                                    "HomDel", "homdel"),collapse = "|"),
    cna_incompletedel_annotation = c("DEL", "del", "Del"),
    chrom_names = c(
      c(seq_len(23), "X", "Y"),
      paste0("chr", c(seq_len(23), "X", "Y"))
    ),
    colnames_gene = c("GENE", "gene", "Gene", 
                      "Gene_name", "GENE_NAME", "gene_name"),
    colnames_af = c("AF", "af", "Af"),
    colnames_ref = c("REF", "ref", "Ref"),
    colnames_alt = c("ALT", "alt","Alt"),
    colnames_tcn = c("TCN", "tcn", "Tcn"),
    colnames_cna_type = c("cna_type", "CNA_type", "Cna_Type", "CNA_Type"),
    colnames_all_imb = c("allelic_imbalance", "all_imb", "AI"),
    colnames_genotype = c("gt_cna", "genotype", "GT", "gt"),
    sex_names = c("male", "m", "female", "f") %>% c(.,toupper(.))
  )
  res <- type_list[[which_one]]
  if(regex){
    res_regex <- paste0("^", res, "$")
    return(res_regex)
  } else {
    return(res)
  }
  
}
