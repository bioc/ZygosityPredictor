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
#' @importFrom IRanges subsetByOverlaps
#' @export
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
                                is_pre_eval=TRUE
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
  geneModel <- check_gr_gene_model(geneModel, is_pre_eval)
  assumeSomCnaGaps <- check_opt_assgap(assumeSomCnaGaps, ploidy)
  somCna <- check_somCna(somCna, geneModel, sex, ploidy, 
                         assumeSomCnaGaps, colnameTcn, 
                         colnameCnaType) 
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
  ## predict zygosity for each gene
  if(!is.null(df_germ)|!is.null(df_som)|!is.null(df_homdels)){
    df_all_mutations <- combine_main_variant_tables(df_germ, df_som, 
                                                    df_homdels,
                                                    templateGenes, purity)
  }
  if(is_pre_eval==TRUE){
    return(list(evaluation_per_variant= 
                  bind_incdel_to_pre_eval(df_incompletedels, df_all_mutations),
                combined_uncovered=combined_uncovered
                
    ))
  } else {
    return(list(evaluation_per_variant=df_all_mutations,
                df_incompletedels=df_incompletedels,
                combined_uncovered=combined_uncovered))    
  }
  
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
                             haploBlocks=NULL,
                             phasedVcf=NULL,
                             distCutOff=5000,
                             verbose=FALSE,
                             logDir=NULL,
                             refGen="hg19"
){
  status <- info <- wt_cp <- . <- df_homdels <- df_all_mutations <- gene <-
    final_phasing_info <- combined_read_details <-  final_output <-
    uncovered_som <- uncovered_germ <- gr_germ_cov <- gr_som_cov <- 
    som_covered <- germ_covered <- final_phasing_info <- 
    combined_snp_phasing <- evaluation_per_gene <- NULL
  #vm("initializing evaluation per variant", verbose)
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
                                                    byTcn,
                                                    is_pre_eval=FALSE)
  evaluation_per_variant <- evaluation_per_variant_pre$evaluation_per_variant
  if(!is.null(evaluation_per_variant)){
    df_all_mutations <- evaluation_per_variant
    if(!nrow(df_all_mutations)==0){
      bamDna <- check_bam(bamDna)
      bamRna <- check_rna(bamRna)
      logDir <- check_logDir(logDir)
      vcf <- check_vcf(vcf) 
      haploBlocks <- check_haploblocks(haploBlocks)
      phasedVcf <- check_vcf(phasedVcf)
      #vm("initializing evaluation per gene", verbose)
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
        haploBlocks,
        phasedVcf,
        distCutOff, 
        verbose,
        logDir,
        refGen)
      #print(result_list)
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
        # combined_read_details <- lapply(result_list, nth, n=3) %>% compact() %>%
        #   bind_rows()
        # combined_snp_phasing <- lapply(result_list, nth, n=4) %>% compact() %>%
        #   bind_rows()
      } 
    }
    evaluation_per_gene <- bind_incdel_to_final_eval(
      evaluation_per_variant_pre$df_incompletedels,
      final_output)
  }
  result_list <- list(
    eval_per_variant=evaluation_per_variant,
    eval_per_gene=evaluation_per_gene,
    phasing_info=final_phasing_info,
    #readpair_info=combined_read_details,
    uncovered_input=evaluation_per_variant_pre$combined_uncovered#,
    #ext_snp_phasing=combined_snp_phasing
  ) %>%
    compact()
  return(result_list)
} 