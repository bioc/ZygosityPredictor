# library(testthat)
# library(dplyr)
# library(stringr)
# library(GenomicRanges)
# 
# setwd("C:/Users/macrh/repos/ZygosityPredictor")
# source("R/fncts.R")
# 
# test_that(
#   "Test ZygosityPredictor functions",
#   {
#     expect_equal(aff_som_copies(3, 0.5, 2, 1, "f"), 1)
#     expect_equal(aff_som_copies("X", 0.5, 2, 0.5, "m"), 1.5)
#     expect_equal(aff_som_copies("X", 0.5, 2, 0.5, "f"), 2)
#     
#     expect_equal(aff_germ_copies("X", 0.5, 2, 0.5, "f"), 1)
#     expect_equal(aff_germ_copies("X", 0.5, 2, 0.5, "m"), 0.5)    
#     
#     expect_equal(check_tcn(2), 2)
#     expect_error(check_tcn("two"))
#     expect_equal(check_chr(22), "22")
#     expect_equal(check_chr("chr22"),"chr22")
#     expect_error(check_chr("chr44"))
#     expect_error(check_somCna("I am a sCNA result"))
#     expect_error(check_gr_gene_model("I am a gene model"))
#     expect_null(check_gr_som(NULL))
#     expect_null(check_gr_germ(NULL))
#     expect_equal(check_purity("0.1"), 0.1)
#     expect_error(check_purity(4))
#     expect_null(check_ploidy(NULL))
#     expect_equal(check_ploidy("10"), 10)
#     expect_equal(check_sex("F"), "female")
#     
#   }
# )
# 
# 
