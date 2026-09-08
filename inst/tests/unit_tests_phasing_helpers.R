library(testthat)
library(ZygosityPredictor)

## ---------------------------------------------------------------------------
## Unit tests for the pure helper functions used throughout the haplotype-
## phasing pipeline (R/phasing.R and R/general.R).
##
## These tests define the "contract" that must be preserved when optimizing
## the hot-path functions. If an optimization changes the observable output
## of any of these functions, the corresponding test here will fail.
## ---------------------------------------------------------------------------

test_that("get_single_index and get_xy_index are exact inverses", {

  nr <- 5
  ## indices in row-major order (column index cycle, row-major storage)
  single <- ZygosityPredictor:::get_single_index(2, 3, nr)  # x=2, y=3
  expect_equal(single, nr * (2 - 1) + 3)  # 5 + 3 = 8

  ## a vector of single indices
  inp <- c(1, 2, 5, 8, 25)
  xy <- ZygosityPredictor:::get_xy_index(inp, nr)

  expect_equal(nrow(xy), length(inp))
  expect_equal(ncol(xy), 2)
  expect_true(all(colnames(xy) %in% c("x", "y")))

  ## round-trip: single -> xy -> single must reproduce the input
  back <- apply(xy, 1, function(r) ZygosityPredictor:::get_single_index(r["x"], r["y"], nr))
  expect_equal(as.numeric(back), inp)

  ## edge case: exact multiple of nr (y wraps to nr)
  expect_equal(as.numeric(ZygosityPredictor:::get_xy_index(5, 5)["y"]), 5)
  expect_equal(as.numeric(ZygosityPredictor:::get_xy_index(5, 5)["x"]), 1)

  ## last element of matrix (x=nr, y=nr)
  expect_equal(as.numeric(ZygosityPredictor:::get_xy_index(25, 5)["x"]), 5)
  expect_equal(as.numeric(ZygosityPredictor:::get_xy_index(25, 5)["y"]), 5)
})

test_that("get_xy_index handles the full index range of a square matrix", {
  nr <- 4
  all_inp <- 1:(nr * nr)
  xy <- ZygosityPredictor:::get_xy_index(all_inp, nr)
  back <- apply(xy, 1, function(r) ZygosityPredictor:::get_single_index(r["x"], r["y"], nr))
  expect_equal(as.numeric(back), all_inp)
})

test_that("get_xy_index returns correct row-major coordinates", {
  nr <- 3
  xy <- ZygosityPredictor:::get_xy_index(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nr)
  ## row-major: first row y=1..3, then x increments
  expect_equal(unname(xy[, "x"]), c(1,1,1,2,2,2,3,3,3))
  expect_equal(unname(xy[, "y"]), c(1,2,3,1,2,3,1,2,3))
})

test_that("make_dist_matrix computes pair distances with distCutOff applied", {
  v <- c(100, 200, 7000, 7500)
  vars <- paste0("m", 1:4)
  mat <- ZygosityPredictor:::make_dist_matrix(v, vars, distCutOff = 5000)

  expect_equal(dim(mat), c(4, 4))
  expect_equal(rownames(mat), vars)
  expect_equal(colnames(mat), vars)

  ## lower triangle set to 0
  expect_true(all(mat[lower.tri(mat)] == 0))

  ## diagonal (self-distance) is 0
  expect_true(all(diag(mat) == 0))

  ## within cutoff: |200-100| = 100 stored in the UPPER triangle
  expect_equal(mat[1, 2], 100)
  ## |7500-7000| = 500, within cutoff -> 500 (upper triangle mat[3,4])
  expect_equal(mat[3, 4], 500)
})

test_that("make_dist_matrix zeros out pairs beyond the cutoff", {
  v <- c(100, 200, 10000)
  mat <- ZygosityPredictor:::make_dist_matrix(v, paste0("m", 1:3), distCutOff = 5000)
  ## |10000-100| = 9900 > 5000 -> 0 (upper triangle mat[1,3])
  expect_equal(mat[1, 3], 0)
  ## |10000-200| = 9800 > 5000 -> 0
  expect_equal(mat[2, 3], 0)
  ## |200-100| = 100 -> kept
  expect_equal(mat[1, 2], 100)
})

test_that("ascii_to_dec converts Phred quality letters to probabilities", {
  ## Phred+33: '!' has ASCII 33 -> Q=0 -> error prob 10^(0/-10)=1
  expect_equal(ZygosityPredictor:::ascii_to_dec("!"), 1)
  ## 'I' has ASCII 73 -> Q=40 -> 10^(-4) = 0.0001
  expect_equal(ZygosityPredictor:::ascii_to_dec("I"), 10^(-40/10))
  ## NA stays NA
  expect_true(is.na(ZygosityPredictor:::ascii_to_dec(NA)))
})

test_that("ascii_to_dec averages over multiple quality letters", {
  ## two identical letters -> mean of identical = same value
  expect_equal(ZygosityPredictor:::ascii_to_dec("II"), 10^(-40/10))
  ## 'I' (Q40 -> 1e-4) and '!' (Q0 -> 1) -> mean of exponents
  ## exp = (Q-33)/-10 -> I: -4, !: 0 -> mean -2 -> 10^-2 = 0.01
  expect_equal(ZygosityPredictor:::ascii_to_dec("I!"), 0.01)
})

test_that("define_class classifies variants correctly", {
  ## SNV: equal length
  expect_equal(ZygosityPredictor:::define_class("A", "G"), "snv")
  ## insertion: alt longer
  expect_equal(ZygosityPredictor:::define_class("A", "AG"), "ins")
  ## deletion: ref longer
  expect_equal(ZygosityPredictor:::define_class("AG", "A"), "del")
})

test_that("get_classification respects homdel precedence", {
  data <- data.frame(class = "homdel", alt = "A", ref = "G")
  expect_equal(ZygosityPredictor:::get_classification(data), "homdel")
  data2 <- data.frame(class = "snv", alt = "A", ref = "G")
  expect_equal(ZygosityPredictor:::get_classification(data2), "snv")
})

test_that("calc_left_wt_copies follows the phasing combination logic", {
  mtcn <- 2
  ## same-phased (nconst==1): left wt = mtcn - max(aff1, aff2)
  expect_equal(ZygosityPredictor:::calc_left_wt_copies(mtcn, 1, 0.5, 0.5), 1.5)
  ## diff-phased (nconst==2): left wt = mtcn - sum(aff1, aff2)
  expect_equal(ZygosityPredictor:::calc_left_wt_copies(mtcn, 2, 0.5, 0.5), 1.0)
  expect_equal(ZygosityPredictor:::calc_left_wt_copies(2, 2, 1.0, 1.0), 0.0)
})

test_that("aggregate_probs combines probabilities cumulatively", {
  expect_equal(ZygosityPredictor:::aggregate_probs(0.5), 0.5)
  ## p1 + (1-p1)*p2 = 0.5 + 0.5*0.5 = 0.75
  expect_equal(ZygosityPredictor:::aggregate_probs(c(0.5, 0.5)), 0.75)
  ## 0.5 + 0.5*0.5 + 0.75*0.5 = 0.5+0.25+0.375 = 1.125... but formula is iterative:
  ## np=0.5; np=0.5+(1-0.5)*0.5=0.75; np=0.75+(1-0.75)*0.5=0.875
  expect_equal(ZygosityPredictor:::aggregate_probs(c(0.5, 0.5, 0.5)), 0.875)
})

test_that("get_string_const maps nconst to const labels", {
  expect_equal(ZygosityPredictor:::get_string_const(2), "diff")
  expect_equal(ZygosityPredictor:::get_string_const(1), "same")
  expect_equal(ZygosityPredictor:::get_string_const(0), "null")
  expect_equal(ZygosityPredictor:::get_string_const(NA), "null")
})

test_that("make_dist_matrix input validation for empty/single vectors", {
  ## single variant: 1x1 matrix of 0
  mat <- ZygosityPredictor:::make_dist_matrix(100, "m1", 5000)
  expect_equal(dim(mat), c(1, 1))
  expect_equal(mat[1, 1], 0)
})

test_that("formula_genotype_likelihood returns expected structure", {
  res <- ZygosityPredictor:::formula_genotype_likelihood(
    gtl_g = 2, gtl_eps_l = c(0.01, 0.01),
    gtl_eps_v = c(0.01, 0.01),
    gtl_m = 2, gtl_k = 4
  )
  expect_equal(length(res), 4)
  expect_named(res, c("gt", "lklhd", "prod_ref", "prod_alt"))
  expect_equal(as.numeric(res["gt"]), 2)
})

test_that("get_main_mut_pos returns main mutation indices", {
  ## set up a minimal ZP_env with main muts and a phasing matrix
  ZP_env <- new.env()
  ZP_env$global_ZygosityPredictor_variable_main_muts <- c("m1", "m2", "m3")
  ZP_env$global_ZygosityPredictor_variable_mat_phased <-
    matrix(0, nrow = 3, ncol = 3)
  ## main mut positions: 1,2,3 -> seq(1,3) for first block etc.
  ## with lm=3, nrow=3: block i covers seq((i-1)*3+1, (i-1)*3+3)
  ## actually get_main_mut_pos returns for each of 1..lm a seq of lm indices
  res <- ZygosityPredictor:::get_main_mut_pos(ZP_env)
  ## for lm=3: lapply over 1:3 -> seq((i-1)*3+1, (i-1)*3+3)
  expect_equal(res, c(1,2,3, 4,5,6, 7,8,9))
})

test_that("get_main_mut_conns counts connections per main mutation", {
  ZP_env <- new.env()
  ZP_env$global_ZygosityPredictor_variable_main_muts <- c("m1", "m2")
  ZP_env$global_ZygosityPredictor_variable_mat_phased <- matrix(0, 2, 2)
  rownames(ZP_env$global_ZygosityPredictor_variable_mat_phased) <- c("m1", "m2")
  colnames(ZP_env$global_ZygosityPredictor_variable_mat_phased) <- c("m1", "m2")
  ## fill with a connection: m1-m2
  ZP_env$global_ZygosityPredictor_variable_mat_phased[1, 2] <- 1
  conns <- ZygosityPredictor:::get_main_mut_conns(ZP_env)
  expect_true(all(c("m1", "m2") %in% names(conns)))
})
