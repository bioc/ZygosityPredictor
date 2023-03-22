library(testthat)
library(ZygosityPredictor)
#library(GenomicRanges)
## prepare inputs for predict_zygosity() unit test

FILE_BAM <- "inst/extdata/ZP_example.bam"#
PURITY <- 0.98
PLOIDY <- 1.57
SEX <- "female"

GR_SCNA  = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(170060, 34520990),
    end = c(34520990, 83198614),
    tcn = c(2, 1),
    cna_type = c("neutral", "LOH")
  )
)
GR_SOM_SMALL_VARS = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 7675088,
    end = 7675088,
    ref = "C",
    alt = "T",
    af = 0.65,
    gene = "TP53"
  )
)
GR_GERM_SMALL_VARS = GenomicRanges::GRanges(
  dplyr::tibble(
    chr="chr17",
    start = 41771694,
    end = 41771694,
    ref = "GTGT",
    alt = "G",
    af = 0.95,
    gene = "JUP"
  )
)
GR_GENE_MODEL = GenomicRanges::GRanges(
  dplyr::tibble(
    chr = "chr17",
    start = c(7661778, 41754603),
    end = c(7687538, 41786931),
    gene = c("TP53", "JUP")
  )
)

# predict_zygosity(purity = purity, sex = sex,
#   somCna = cnvs,
#   somSmallVars = somatic_vars,
#   germSmallVars = germline_vars,
#   geneModel = reference,
#   bamDna = bamfile
# )


test_that(
  "Test ZygosityPredictor exportet functions",
  {
    expect_equal(aff_som_copies(3, 0.5, 2, 1, "f"), 1)
    expect_equal(aff_som_copies("X", 0.5, 2, 0.5, "m"), 1.5)
    expect_equal(aff_som_copies("X", 0.5, 2, 0.5, "f"), 2)

    expect_equal(aff_germ_copies("X", 0.5, 2, 0.5, "f"), 1)
    expect_equal(aff_germ_copies("X", 0.5, 2, 0.5, "m"), 0.5)

    expect_type(predict_zygosity(
      purity = PURITY,
      ploidy = PLOIDY,
      sex = SEX,
      somCna = GR_SCNA,
      somSmallVars = GR_SOM_SMALL_VARS,
      germSmallVars = GR_GERM_SMALL_VARS,
      geneModel = GR_GENE_MODEL,
      bamDna = FILE_BAM
    ), "list")

    expect_message(predict_zygosity(
      purity = PURITY,
      ploidy = PLOIDY,
      sex = SEX,
      somCna = GR_SCNA,
      geneModel = GR_GENE_MODEL,
      bamDna = FILE_BAM
    )
    )

    expect_error(predict_zygosity(
      purity = PURITY,
      ploidy = PLOIDY,
      sex = SEX,
      somCna = GR_SCNA[,1],
      somSmallVars = GR_SOM_SMALL_VARS,
      germSmallVars = GR_GERM_SMALL_VARS,
      geneModel = GR_GENE_MODEL,
      bamDna = FILE_BAM
    ))
  }
)

test_that(
  "Test ZygosityPredictor internal functions",
  {
    expect_equal(ZygosityPredictor:::check_tcn(2), 2)
    expect_error(ZygosityPredictor:::check_tcn("two"))
    expect_equal(ZygosityPredictor:::check_chr(22), "22")
    expect_equal(ZygosityPredictor:::check_chr("chr22"),"chr22")
    expect_error(ZygosityPredictor:::check_chr("chr44"))
    expect_error(ZygosityPredictor:::check_somCna("I am a sCNA result"))
    expect_error(ZygosityPredictor:::check_gr_gene_model("I am a gene model"))
    expect_null(ZygosityPredictor:::check_gr_small_vars(NULL, origin="somatic"))
    expect_equal(ZygosityPredictor:::check_purity("0.1"), 0.1)
    expect_error(ZygosityPredictor:::check_purity(4))
    expect_null(ZygosityPredictor:::check_ploidy(NULL))
    expect_equal(ZygosityPredictor:::check_ploidy("10"), 10)
    expect_equal(ZygosityPredictor:::check_sex("F"), "female")

  }
)


