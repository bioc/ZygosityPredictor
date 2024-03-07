

library(tidyverse)
library(MultiAssayExperiment)
library(RaggedExperiment)

#c("H021-1XUB9H_tumor_control02_WGS", "H021-25LC94_tumor06_blood02_WGS", "H021-26BUQE_tumor_blood_WGS")
load_ZP_functions <- function(){
  source("/home/m168r/projects/ZygosityPredictor/R/check_input.R")
  source("/home/m168r/projects/ZygosityPredictor/R/exported.R")
  source("/home/m168r/projects/ZygosityPredictor/R/general.R")
  source("/home/m168r/projects/ZygosityPredictor/R/phasing.R")
}

file_paths <- tibble(
  ## reduced dataMASTER object to 5 test cases
  red_dataMASTER = "/home/m168r/home_extension/test_tools/dataMASTER_reduced_env",
  ## path to dataMASTER object
  dataMASTER = "/omics/odcf/analysis/OE0246_projects/datamaster/dataMASTER_release/object/dataMASTER.RData",
  ## path to TOP-ART calculator helper functions... cloned repo from gitlab
  TAC_preprocessing = "/home/m168r/projects/TopArtCalc/TAC_preprocessing/pre_proc_fncts.R",
  ## path to TAC file selection scipts
  TAC_file_selection="/home/m168r/projects/TopArtCalc/TAC_file_selection/file_sel_fncts.R",
  ## path to hg19 gene model
  file_gene_model = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_Exons_plain.bed.gz",
)

load(file_paths$red_dataMASTER)

## load gene model
GR_GENE_MODEL <- read_tsv(
  file_paths$file_gene_model,
  col_names=c("chr", "start", "end","gene", "n_exon", "strand", "cdsstart", "cdsend"),
  skip=1) %>%
  GRanges()


samples_to_run <- 
  c("WGS.1XUB9H.tumor", "WES.25LC94.tumor06", "WGS.26BUQE.tumor", 
    "WES.2WUY1M.tumor02", "WES.2X77QL.metastasis")

source(file_paths$TAC_preprocessing)
source(file_paths$TAC_file_selection)


## genes of interest
goi <- TopArtCalc::load_somatic_topart_genes()

## the user can control the genes to be assessed by ZygosiytPredictor by filtering
## the gene model:
geneModel <- GR_GENE_MODEL[which(GR_GENE_MODEL$gene %in% goi)]


#library(ZygosityPredictor)



SAMPLE_ID <- "WGS.1XUB9H.tumor"
myobj = new_obj[,SAMPLE_ID,]

prepared_input <- prepare_ZP_dataMASTER(SAMPLE_ID, myobj)


load_ZP_functions()


full_prediction <- predict_zygosity(
  purity=colData(myobj)$TumorCellContent,
  sex=colData(myobj)$Sex,
  somCna = prepared_input$gr_cna,
  geneModel = geneModel,
  bamDna = prepared_input$bamDna,
  somSmallVars = prepared_input$gr_som,
  germSmallVars = prepared_input$gr_germ,
  bamRna = prepared_input$bamRna,
  ploidy = colData(myobj)$Ploidy,
  haploBlocks=prepared_input$haploBlocks,
  vcf=prepared_input$phasedVcf,
  ## options for run
  printLog = F,
  includeIncompleteDel = F,
  includeHomoDel = F,
  verbose=T,
  debug=F,
  ## highest distance of two variants that should be phased (if reduced runtime reduces as well)
  distCutOff = 5000,
  ## if provided, detailed information will be stored there
  logDir=NULL,
  showReadDetail = F,
  ## Secondary phasing approach if read-level phasing fails... can lead to wroing results
  AllelicImbalancePhasing=F
)








ZP_dataMASTER <- lapply(samples_to_run, function(SAMPLE_ID){
  
  print(SAMPLE_ID)
  
  myobj = new_obj[,SAMPLE_ID, ]
  
  prepared_input <- prepare_ZP_dataMASTER(SAMPLE_ID, myobj)
  
  full_prediction <- predict_zygosity(
    purity=colData(myobj)$TumorCellContent,
    sex=colData(myobj)$Sex,
    somCna = prepared_input$gr_cna,
    geneModel = geneModel,
    bamDna = prepared_input$bamDna,
    somSmallVars = prepared_input$gr_som,
    germSmallVars = prepared_input$gr_germ,
    bamRna = prepared_input$bamRna,
    ploidy = colData(myobj)$Ploidy,
    haploBlocks=prepared_input$haploBlocks,
    vcf=prepared_input$phasedVcf,
    ## options for run
    printLog = T,
    includeIncompleteDel = F,
    includeHomoDel = F,
    verbose=F,
    debug=F,
    ## highest distance of two variants that should be phased (if reduced runtime reduces as well)
    distCutOff = 5000,
    ## if provided, detailed information will be stored there
    logDir=NULL,
    showReadDetail = F,
    ## Secondary phasing approach if read-level phasing fails... can lead to wroing results
    AllelicImbalancePhasing=F
  )
  
})
















