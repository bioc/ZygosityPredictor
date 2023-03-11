library(testthat)
library(ZygosityPredictor)

## prepare inputs for predict_zygosity() unit test

FILE_BAM <- "inst/extdata/ZP_example.bam"#
PURITY <- 0.98
PLOIDY <- 1.57
SEX <- "female"
data("GR_SCNA")
data("GR_GERM_SMALL_VARS")
data("GR_SOM_SMALL_VARS")
data("GR_GENE_MODEL")

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


