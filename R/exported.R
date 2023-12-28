#' accesor for gene predictions printing detailed info about how a gene status 
#' was assigned
#' 
#' @param fp full prediction (output of predict_zygoisty())
#' @param inp_gene name of gene that should be printed with detailed information
#' @param n max number of rows to print, as some gene status depend on loads of
#' phasing results#' 
#' 
#' @export
#' @importFrom rlang ensym
gene_ov <- function(fp, inp_gene, n=20){
  #inp_gene <- as.character(rlang::ensym(inp_gene))
  inp_gene <- as.character(ensym(inp_gene))
  eval_per_gene <- fp$eval_per_gene %>% filter(gene==inp_gene) 
  eval_per_var <- fp$eval_per_variant %>% filter(gene==inp_gene) %>%
    mutate_at(.vars=c("af", "tcn", "wt_cp", "aff_cp"),
              .funs=round, 2) %>%
    select(-vn_status, -pre_info, -purity)

  message("Top level: Gene status\n")
  message(print_tibble(eval_per_gene))
  message("\n\nSub level: Evaluation per variant\n")
  message(print_tibble(eval_per_var))
  
  if(!is.null(fp$phasing_info)){
   
    phasing_info <- fp$phasing_info %>% 
      filter(gene==inp_gene) %>%
      mutate_at(.vars=c("conf", "wt_cp", "min_poss_wt_cp", "max_poss_wt_cp"),
                .funs=round, 2) %>% 
      arrange(desc(nconst), desc(conf))
    if(!is.null(fp$detailed_RLP_info)){
      detailed_RLP <- fp$detailed_RLP_info %>% 
        filter(gene==inp_gene) %>%
        select(-conf, -xsq_diff, -xsq_same, -v_same, -v_diff, -phasing) %>%
        mutate_at(.vars=c("mut2", "mut1", "both", "p_same",    "p_diff"),
                  .funs=round, 2) %>%
        arrange(desc(nconst))   
      if(nrow(detailed_RLP)<n){
        n_RLP <- nrow(detailed_RLP)
      } else {
        n_RLP <- n
      }  
      RLP_mes <- paste(c("\n\nSub level: All read-level-phasing combinations, including SNPs\n", "Showing ", n_RLP, " of ", nrow(detailed_RLP), " phasing attempts\n\n",
          print_tibble(detailed_RLP[1:n_RLP,]), "\n"), collapse="")
    } else {
      RLP_mes <- ""
    }
    if(!is.null(fp$detailed_AIP_info)){
      detailed_AIP <- fp$detailed_AIP_info %>% 
        filter(gene==inp_gene) %>%
        #select(-conf, -xsq_diff, -xsq_same, -v_same, -v_diff, -phasing) %>%
        # mutate_at(.vars=c("p_same",    "p_diff"),
        #           .funs=round, 2) %>%
        arrange(desc(nconst))
  
      if(nrow(detailed_AIP)<n){
        n_AIP <- nrow(detailed_AIP)
      } else {
        n_AIP <- n
      }
      AIP_mes <- paste(c("\n\nSub level: All allelic-imbalance-phasing combinations, including SNPs\n", "Showing ", n_AIP, " of ", nrow(detailed_RLP), " phasing attempts\n\n",
                         print_tibble(detailed_RLP[1:n_AIP,]), "\n"), collapse="")
    } else {
      AIP_mes <- ""
    }
    if(nrow(phasing_info)>0){
      message("\n\nSub level: Main phasing combinations\n")
      message(print_tibble(phasing_info))
      message(RLP_mes)
      message(AIP_mes)
      # message("\n\nSub level: All read-level-phasing combinations, including SNPs\n", "Showing ", n, " of ", nrow(detailed_RLP), " phasing attempts\n")
      # message(print_tibble(detailed_RLP[1:n,]))
      # message("\n\nSub level: All allelic-imbalance-phasing combinations, including SNPs\n", "Showing ", n, " of ", nrow(detailed_AIP), " phasing attempts\n")
      # message(print_tibble(detailed_AIP[1:n,]))
    }
  }

  
}
#' accesor for ZygoistyPredictor runs. Prints an overview about the run
#' 
#' @param fp full prediction (output of predict_zygoisty())
#' 
#' @export
ZP_ov <- function(fp){
  
  if(!is.null(fp$eval_per_variant)){
    n_suc_vars <- nrow(fp$eval_per_variant)
    if(!is.null(fp$uncovered_input)){
      n_uncovered_vars <- nrow(fp$uncovered_input)
    } else {
      n_uncovered_vars <- 0
    }
    n_vars <- n_suc_vars+n_uncovered_vars
    if(!is.null(fp$eval_per_gene)){
      n_pred_genes <- nrow(fp$eval_per_gene)
      overview_pred_genes <- fp$eval_per_gene %>%
        group_by(status, eval_by, .drop=FALSE) %>%
        tally() %>%
        print_tibble()
    } else {
      n_pred_genes <- 0
    }
    if(!is.null(fp$phasing_info)){
      phased_genes <- fp$phasing_info %>%
        filter(nconst>0) %>% pull(gene) %>%
      
        unique() %>% sort %>% paste(collapse = ", ")
    }
    message("Zygosity Prediction of ", n_vars, " input variants")
    message(n_uncovered_vars, " are uncovered by sCNA input and were not evaluated")
    message(n_pred_genes, " genes analyzed:")
    message(overview_pred_genes)
    message("phased genes:\n", phased_genes)
  } else {
    message("No input variants provided to ZygosityPredictor")
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
#' @param af_normal Allele-frequency in normal tissue (numeric value between 0 
#' and 1) 0.5 represents heterozygous variants in diploid genome, 
#' 1 would be homozygous. Could be relevant if germline CNVs are present at the 
#' position. Then also the c_normal parameter would have to be adjusted.
#' @return A numeric value indicating the affecting copies for the variant
#' @export
#' @examples 
#' library(dplyr)
#' library(purrr)
#' library(stringr)
#' aff_germ_copies(af=0.67, tcn=2, purity=0.9, chr="chrX", sex="female")
aff_germ_copies <- function(chr, af, tcn, purity, sex, 
                            c_normal=NULL, af_normal=0.5){
  cv <- formula_checks(chr, af, tcn, purity, sex, c_normal, af_normal)
  aff_cp <- cv$af*cv$tcn+(cv$af-cv$af_normal)*cv$c_normal*((1/cv$purity)-1)
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
#' predicts zygosity of a set of variants
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
#' @param includeHomoDel default = TRUE; if FALSE homozygous deleteions are 
#' excluded
#' @param includeIncompleteDel default = TRUE; if FALSE heterzygous deleteions 
#' are excluded
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
#' @param verbose logical, default=FALSE; prints functions that are called
#' @importFrom IRanges subsetByOverlaps
#' @export
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
#' predict_per_variant(purity = purity, sex = sex, 
#'   somCna = cnvs,
#'   somSmallVars = somatic_vars,
#'   germSmallVars = germline_vars,
#'   geneModel = reference
#' )
predict_per_variant <- function(purity, 
                                sex,
                                somCna, 
                                geneModel=NULL,
                                somSmallVars=NULL, 
                                germSmallVars=NULL,
                                ploidy=NULL, 
                                colnameTcn=NULL,
                                colnameCnaType=NULL,
                                includeHomoDel=TRUE,
                                includeIncompleteDel=TRUE,
                                assumeSomCnaGaps=FALSE,
                                byTcn=TRUE,
                                verbose=FALSE
){
  status <- info <- wt_cp <- . <- df_homdels <- df_all_mutations <- gene <-
    final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- 
    combined_snp_phasing <- NULL
  if(!exists("global_ZygosityPredictor_variable_embedded")){
    is_pre_eval <- TRUE
    set_global_variables(FALSE, verbose, FALSE)
    somCna <- check_somCna(somCna, geneModel, sex, ploidy, 
                           assumeSomCnaGaps, colnameTcn, 
                           colnameCnaType)  
    purity <- check_purity(purity)
    sex <- check_sex(sex)
  } else {
    is_pre_eval <- FALSE
  }
  func_start()
  ## check input for valid format
  somSmallVars <- check_gr_small_vars(somSmallVars, "somatic")
  germSmallVars <- check_gr_small_vars(germSmallVars, "germline")
  ploidy <- check_ploidy(ploidy)
  geneModel <- check_gr_gene_model(geneModel, is_pre_eval)
  assumeSomCnaGaps <- check_opt_assgap(assumeSomCnaGaps, ploidy)
  includeIncompleteDel <- check_opt_incdel(includeIncompleteDel, ploidy)
  includeHomoDel <- check_opt_incdel(includeHomoDel, ploidy)
  if(is.null(geneModel)){
    ## only small variants can be evaulated
    if(includeIncompleteDel==TRUE|includeHomoDel==TRUE){
      includeIncompleteDel <- FALSE
      includeHomoDel <- FALSE
      warning("To include large deletions the geneModel muist be provided.",
              "IncludeHomoDel and IncludeIncompleteDel will be FALSE")
    }
    templateGenes <- c(somSmallVars$gene, germSmallVars$gene)
  } else {
    templateGenes <- geneModel$gene
  }
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
  df_incompletedels <- extract_all_dels_of_sample(somCna, geneModel, 
                                                  "incompletedel", TRUE, sex, 
                                                  ploidy, includeIncompleteDel)
  
  #print(df_germ)
  #print(df_som)
  #print()
  
  ## predict zygosity for each gene
  if(!is.null(df_germ)|!is.null(df_som)|!is.null(df_homdels)){
    #print(1)
    df_all_mutations <- combine_main_variant_tables(df_germ, df_som, 
                                                    df_homdels,
                                                    templateGenes, purity)
    #print(2)
  }
  if(is_pre_eval){
    
    full_eval <- list(evaluation_per_variant= 
                  bind_incdel_to_pre_eval(df_incompletedels, df_all_mutations),
                combined_uncovered=combined_uncovered)
  } else {
    full_eval <- list(evaluation_per_variant=df_all_mutations,
                df_incompletedels=df_incompletedels,
                combined_uncovered=combined_uncovered)   
  }
  func_end()
  if(is_pre_eval){
    remove_global_vars()
  }
  return(full_eval)
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
#' @param AllelicImbalancePhasing logical, default=FALSE. Enables alleleic imbalance phasing if TRUE
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
#' @param debug logical, default=FALSE; prints output for debugging
#' @param verbose logical, default=FALSE; prints functions that are called
#' @param haploBlocks GRanges object containing haploblocks. Haploblocks are
#' defined as genomic regions in which SNPs are phased to a specific allele.
#' For example a haploblock could be chr1:1000-10000. This would mean that every
#' genotype annotation in the format "1|0" or "0|1" of a SNP in this region will 
#' be used to phase somatic variants and define their genotype
#' @param logDir character; path to directory where logfiles and detailed infos 
#' of the run can be stored, if not given, no details will be stored or printed     
#' @param snpQualityCutOff numeric, default=1; Cutoff to filter for SNPS that 
#' can be used for phasing                     
#' @param phasingMode character, default="fast"; if set to full. Even if high 
#' confidence phasing result could be achieved, following phasing approaches 
#' will be carried out
#' @param showReadDetail default = FALSE; if TRUE a table is added to the 
#' output, containing all used reads/rea-pairs with anntated read classification 
#' (mut1, mut2, both, none, skipped, dev_var)
#' 
#' 
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
                             haploBlocks=NULL,
                             distCutOff=5000,
                             verbose=FALSE,
                             debug=FALSE,
                             logDir=NULL,
                             snpQualityCutOff=1, 
                             phasingMode="fast",
                             AllelicImbalancePhasing=FALSE
){
  status <- info <- wt_cp <- . <- df_homdels <- evaluation_per_variant <- 
    gene <- final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- detailed_RLP_info <- 
    detailed_AIP_info <- 
    combined_snp_phasing <- evaluation_per_gene <- log_list_per_gene <- 
    detailed_phasing_info <- comb_mat_phased <- comb_mat_info <- 
    combined_eval_per_gene <- full_phasing_info <- NULL
  ## define global debugging variable
  set_global_variables(debug, verbose, printLog)
  global_ZygosityPredictor_variable_embedded <<- TRUE
  func_start()
  somCna <- check_somCna(somCna, geneModel, sex, ploidy, 
                          assumeSomCnaGaps, colnameTcn, 
                          colnameCnaType)
  purity <- check_purity(purity)
  sex <- check_sex(sex)
  evaluation_per_variant_pre <- predict_per_variant(purity, 
                                                    sex,
                                                    somCna, 
                                                    geneModel,
                                                    somSmallVars, 
                                                    germSmallVars,
                                                    ploidy, 
                                                    colnameTcn,
                                                    colnameCnaType,
                                                    includeHomoDel,
                                                    includeIncompleteDel,
                                                    assumeSomCnaGaps,
                                                    byTcn)
  evaluation_per_variant <- evaluation_per_variant_pre$evaluation_per_variant
  if(!is.null(evaluation_per_variant)){
    if(!nrow(evaluation_per_variant)==0){
      bamDna <- check_bam(bamDna)
      bamRna <- check_rna(bamRna)
      logDir <- check_logDir(logDir)
      haploBlocks <- check_haploblocks(haploBlocks)
      vcf <- check_vcf(vcf)
      per_gene <- lapply(
        unique(evaluation_per_variant$gene), 
        predict_zygosity_genewise, 
        evaluation_per_variant, 
        bamDna,
        bamRna,
        showReadDetail,
        printLog,
        purity,
        sex,
        haploBlocks,
        vcf,
        distCutOff, 
        logDir,
        somCna,
        snpQualityCutOff, 
        phasingMode,
        AllelicImbalancePhasing)
      full_eval_per_gene <- lapply(per_gene, nth, n=2) %>% compact()
      log_list_per_gene <- lapply(per_gene, nth, n=3)
      if(length(full_eval_per_gene)!=0){
        full_phasing_info <- lapply(full_eval_per_gene, nth, n=2) %>% 
          compact() %>% 
          bind_rows() 
        detailed_RLP_info <- lapply(full_eval_per_gene, nth, n=3) %>% 
          compact() %>% 
          bind_rows() 
        detailed_AIP_info <- lapply(full_eval_per_gene, nth, n=6) %>% 
          compact() %>% 
          bind_rows() 
        comb_mat_phased <- lapply(full_eval_per_gene, nth, n=4) %>%
          compact() %>%
          Reduce(function(x,y)append(x,y),.)
        comb_mat_info <- lapply(full_eval_per_gene, nth, n=5) %>%
          compact() %>%
          Reduce(function(x,y)append(x,y),.)
        combined_eval_per_gene <- lapply(full_eval_per_gene, nth, n=1) %>% 
          bind_rows() %>%
          select(gene, n_mut, status, conf, eval_by, wt_cp, 
                 warning, wt_cp_range, info,  phasing, eval_time_s)
      } 
    }
    evaluation_per_gene <- bind_incdel_to_final_eval(
      evaluation_per_variant_pre$df_incompletedels,
      combined_eval_per_gene) %>%
      transform(status=factor(status, 
                              levels=c("all_copies_affected", "wt_copies_left", "undefined")))
    
  }
  result_list <- list(
    eval_per_variant=evaluation_per_variant,
    eval_per_gene=evaluation_per_gene,
    phasing_info=full_phasing_info,
    uncovered_input=evaluation_per_variant_pre$combined_uncovered,
    log_list_per_gene=log_list_per_gene,
    detailed_RLP_info=detailed_RLP_info,
    detailed_AIP_info=detailed_AIP_info,
    mat_phased=comb_mat_phased,
    mat_info=comb_mat_info
  ) %>%
    compact()
  remove_global_vars()
  func_end()
  return(result_list)
} 