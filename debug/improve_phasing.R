library(igraph)
library(GenomicAlignments)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)
library(vcfR)
library(knitr)

#library(ZygosityPredictor)

## check which samples have genes that cannot be phased

#all_phasing_files <- list.files("/omics/odcf/analysis/hipo/hipo_021/cohort_analysis/TOP-ART/results_per_pid/", pattern="main_phasing", recursive = T, full.names = T)

#loaded <- lapply(all_phasing_files, read_tsv) %>% bind_rows()
#read_tsv(all_phasing_files[1])
#write_tsv(loaded, file=file.path("/home/m168r/home_extension/ZP_revision/phasing_results.tsv"))

#fp <- read_tsv("/home/m168r/home_extension/ZP_revision/phasing_results.tsv")

#poss <- loaded %>%
 # filter(status=="null",
  #       dist<1000)


source("/home/m168r/projects/TopArtCalc/workflow/runTAC_fncts.R")

sf <- read_tsv("/omics/odcf/analysis/hipo/hipo_021/cohort_analysis/TOP-ART/intermediate_data/MASTER/selected_files_PIDSAM_2022-06-17.tsv") %>%
  mutate(sample_id=paste(PID, sample, seq_method, sep="_"))

#poi <- "H021-2WUY1M_tumor_buffy_coat_WES" ## MUS81 imbalance
#poi <- "H021-GW1GNF_tumor_control_WGS"
#poi <- "H021-3GNGK9_metastasis02_blood_WGS" ## PALB2
poi <- "H021-JRJTPJ_tumor_buffy_coat_WGS"  ## RECQL4
#poi <- "H021-1EEW3X_metastasis_blood_WES"
#poi <- "H021-1HBRFX_metastasis_blood_WES"
#poi <- "H021-QGQNDG_metastasis02_buffy_coat_WES"
#poi <- "H021-NXG3Y8_tumor_buffy_coat02_WGS"
#poi <- "H021-D28K7A_metastasis_blood_WGS" ## PRKDC
#poi <- "H021-BDGX7Y_metastasis_blood_WGS" ## PRKDC
#poi <- "H021-SDXUG8_metastasis_blood_WGS" ## ATR
#poi <- "H021-RKAVFH_metastasis_blood_WGS" ## SMG1
#poi <- "H021-Q27DJA_metastasis_buffy_coat02_WGS" ## BRIP1
#poi <- "H021-Y78RB8_metastasis07_blood02_WGS" ## REV3L
#poi <- "H021-ZW1J5Y_metastasis_blood_WGS"
#poi <- "H021-WV44XG_metastasis02_buffy_coat_WES"
#poi <- "H021-VB46RR_tumor_buffy_coat02_WGS" ## POLN
#poi <- "H021-RF6KBF_metastasis_blood_WGS" ## BRCA1
#poi <- "H021-1MDC9C_tumor_blood_WES" ## BRCA1 ext snp phasing
#poi <- "H021-GZV3MZ_metastasis_buffy_coat02_WGS" ## POLM
#poi <- "H021-9QZTUK_tumor_buffy_coat02_WGS"
#poi <- "H021-H1DNZG_tumor_buffy_coat03_WGS"
#poi <- "H021-JRJTPJ_tumor09_buffy_coat_WGS"
#poi <- "H021-HCWEB6_tumor05_blood_WES"
#poi <- "H021-SJEUP3_metastasis_blood_WGS"
#poi <- "H021-FRJ1TB_tumor_buffy_coat04_WGS"
#poi <- "H021-UA82VZ_tumor_control_WGS"
#poi <- "H021-26BUQE_tumor_blood_WGS"
#poi <- "H021-VFUGZW_undefined-neoplasia_blood_WES"
#poi <- "H021-ZSW8VM_metastasis04_buffy_coat02_WGS"
#poi <- "H021-ZR1RTB_tumor_buffy_coat02_WGS"
#poi <- "H021-WK2X94_metastasis03_buffy_coat02_WGS"
#poi <- "H021-GHQ344_tumor02_buffy_coat02_WGS"
#poi <- "P021-AP7A_tumor_blood_WES"#
#poi <- "H021-YFVHVB_metastasis02_blood_WGS"
#poi <- "H021-LDDQVD_tumor_blood_WGS"
#poi <- "H021-Q27DJA_metastasis_buffy_coat02_WGS" ## BRIP2
goi <- c("TIMELESS")
#goi <- TopArtCalc::load_somatic_topart_genes()
files <- sf %>% filter(sample_id == poi)

PID =files$PID
SAMPLE = files$sample
SEQ_METHOD=files$seq_method
cnv_file=files$G1_cnv_file
bamRna=files$rna_file
somatic_snv=files$G1_snv_file
somatic_indel=files$G1_indel_file
charger_file=files$G2_file
bamDna=files$bam_file
FILE_FULL_SNV <- find_raw_snv(files$PID, files$sample, files$seq_method) %>% 
  sort() %>% last()
FILE_FULL_INDEL <- find_raw_indel(files$PID, files$sample, files$seq_method)

phasing_dir <- get_phasing_dir(SAMPLE, SEQ_METHOD, FILE_FULL_SNV)
haploBlocks <- get_haploblock_file(phasing_dir)
phasedVcf <- get_phased_vcfs(phasing_dir)


print("gene_model loaded")

if(SEQ_METHOD=="WES"){
  phasedVcf <- FILE_FULL_SNV
}
vcf <- c(FILE_FULL_SNV, FILE_FULL_INDEL)




SAMPLE_ID <- paste(PID, SAMPLE, SEQ_METHOD, sep="_")

purity <- cnv_file %>% 
  str_match('comb_pro_extra(\\d\\.*\\d*)_(\\d\\.*\\d*).txt$') %>% 
  .[3] %>% as.numeric() 
ploidy <- cnv_file %>% 
  str_match('comb_pro_extra(\\d\\.*\\d*)_(\\d\\.*\\d*).txt$') %>% 
  .[2] %>% as.numeric() 
sex <- get_sex(PID, cnv_file, SEQ_METHOD)
if(!is.null(bamRna)){
  if(is.na(bamRna)|bamRna=="NA"){
    bamRna <- NULL
  } 
}


somCna <- create_cnv_input(cnv_file) 

#goi <- c("POLM")
#goi <- TopArtCalc::load_somatic_topart_genes()
TBL_GENCODE_EXON <- read_tsv("/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_Exons_plain.bed.gz",
                             # n,
                             col_names=c("chr", "start", "end","gene", "exon_number", "strand", "width"),
                             skip=1) %>% #.[19268:nrow(.),]
  filter(gene %in% goi)

somSmallVars <- create_somatic_input(somatic_snv, somatic_indel, PID, TBL_GENCODE_EXON$gene)


geneModel=GenomicRanges::GRanges(TBL_GENCODE_EXON)

germSmallVars <- create_germline_input_new(charger_file, somatic_snv, somatic_indel, TBL_GENCODE_EXON$gene) 
  
if(length(germSmallVars)==0){
  germSmallVars <- NULL
}
  
timestamp <- Sys.time() %>% as.character() %>% str_replace_all(" ", "_") %>% str_replace_all("\\.\\d*$", "")
source("/home/m168r/projects/ZygosityPredictor/R/check_input.R") 
source("/home/m168r/projects/ZygosityPredictor/R/phasing.R")
source("/home/m168r/projects/ZygosityPredictor/R/exported.R")
source("/home/m168r/projects/ZygosityPredictor/R/general.R")  


logDir <- file.path("/home/m168r/home_extension/ZP_revision/test_runs",
                    paste(SAMPLE_ID,
                          timestamp, sep="_"))
unlink(logDir, recursive=T)  
dir.create(logDir)

fp <- predict_zygosity(
    purity=purity,
    sex=sex,
    somCna = somCna,
    geneModel = geneModel,
    bamDna = bamDna,
    somSmallVars = somSmallVars,
    germSmallVars = germSmallVars,
    bamRna = bamRna,
    ploidy = ploidy,
    haploBlocks=haploBlocks,
    vcf=phasedVcf,
    printLog = F,
    includeIncompleteDel = F,
    includeHomoDel = F,
    verbose=F,
    debug=F,
    distCutOff = 5000,
    logDir=logDir,
    showReadDetail = FALSE,
    AllelicImbalancePhasing=T
  )


#################### test quality call
roxygen2::roxygenize("/home/m168r/projects/ZygosityPredictor")
devtools::install_local("/home/m168r/projects/ZygosityPredictor", force=T)

BiocCheck::BiocCheck("/home/m168r/projects/ZygosityPredictor_1.3.2.tar.gz")

colnameTcn=NULL
colnameCnaType=NULL
includeHomoDel=TRUE
includeIncompleteDel=TRUE
showReadDetail=FALSE
printLog=FALSE
assumeSomCnaGaps=FALSE
byTcn=TRUE
distCutOff=5000
verbose=FALSE
refGen="hg19"
snpQualityCutOff <- 1
phasingMode <- "full"
copyNumberPhasing <- T
GENE <- goi

################### load VCF file correctly:





roxygen2::roxygenize("/home/m168r/projects/TopArtCalc/R_tool/")
devtools::install_local("/home/m168r/projects/TopArtCalc/R_tool/", force=T)




cigar <- "11M840048N67M105967N23M"
seq <- "ACAGGTGTCCTGATATGCCTTGTTGCCATGGGATACCTGTTCATGTGTTTTGGAGGCACCGTCTCTCCCTGGGACCAGCTATCGTTCTTCCTCTTCATCAT"
qual <- "AAFF7<JFA<FJJ<<FJFF-<JJ<7FFF-FFJFJ-<J<JA-7A<A7AFAJJ--AJAF-AJJ7F<JFAJFJJ-<7<J-A-A--<JA<<7A<A7F7-F<AAFA"
read_start <-  6574758






################ this for the decison for status
af <- list.files("/home/m168r/home_extension/test_tools/runs/wfZP_latest/results_per_pid/",
                 pattern="direct_phasi", full.names = T, recursive=T) %>%
  lapply(read_tsv, col_types=cols(.default="c")) %>%
  bind_rows()


test <- af %>%
  filter(status!="null") %>%
  mutate_at(.vars = c("sim_diff_x", "sim_diff_p", "sim_same_x", "sim_same_p", "both", "mut1", "mut2"),
            .funs = as.numeric)%>%
  mutate_at(.vars = c("sim_diff_x", "sim_diff_p", "sim_same_x", "sim_same_p"),
            .funs = round, digits=2) %>%
  mutate(status2=case_when(
    sim_diff_x<sim_same_x ~ "diff",
    TRUE ~ "same"
  ),
  dev=status==status2
  ) %>%
  select(1:4, 7, status2, dev, diffx=sim_diff_x, diffp=sim_diff_p, samex=sim_same_x, samep=sim_same_p) %>%
  mutate(
    sum_rel=both+mut1+mut2,
    diffxn=diffx/sum_rel,
    samexn=samex/sum_rel,  
    ratio_x=diffxn/samexn,
    difr_x=abs(diffxn-samexn)
  )









## snv + snv
bamDna <- "/omics/odcf/project/hipo/hipo_021/sequencing/exon_sequencing/view-by-pid/H021-1EEW3X/metastasis/paired/merged-alignment/metastasis_H021-1EEW3X_merged.mdup.bam"



## del + ins
poi <- "H021-1EEW3X_metastasis_blood_WES"
bamDna <- "/omics/odcf/project/hipo/hipo_021/sequencing/exon_sequencing/view-by-pid/H021-1EEW3X/metastasis/paired/merged-alignment/metastasis_H021-1EEW3X_merged.mdup.bam"
ref_chr1 <- "7"
ref_pos1 <- 4821377
ref_gr1 <- GRanges(seqnames = ref_chr1, 
                   ranges = ref_pos1)
ref_chr2 <- "7"
ref_pos2 <- 4821383
ref_gr2 <- GRanges(seqnames = ref_chr2, 
                   ranges = ref_pos2)
ref_alt1 <- "TG"
ref_alt2 <- "C"
ref_ref1 <- "T"
ref_ref2 <- "CAGGT"
ref_class1 <- "ins"
ref_class2 <- "del"


qname <- "ST-K00246:434:HMJY2BBXY:5:1205:14103:29518"





bamDna <- "/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-JRJTPJ/tumor09/paired/merged-alignment/tumor09_H021-JRJTPJ_merged.mdup.bam"
ref_chr1 <- "13"
ref_pos1 <- 108861039
ref_gr1 <- GRanges(seqnames = ref_chr1, 
                   ranges = ref_pos1)

ref_chr2 <- "13"
ref_pos2 <- 108861219
ref_gr2 <- GRanges(seqnames = ref_chr2, 
                   ranges = ref_pos2)


ref_alt1 <- "T"
ref_alt2 <- "T"
ref_ref1 <- "C"
ref_ref2 <- "C"
ref_class1 <- "snv"
ref_class2 <- "snv"

## now load all reads/read-pairs that cover the position of the first variant
#vm("loading reads", verbose)
all_covering_read_pairs <- readGAlignmentPairs(
  bamDna,
  param=ScanBamParam(
    which=GRanges(seqnames = ref_chr1, 
                  ranges = ref_pos1),
    what=c("qname","seq", "cigar", "mapq", "qual")
  )) 


all_reads <- c(
  GenomicAlignments::first(all_covering_read_pairs) %>%
    GRanges(),
  GenomicAlignments::last(all_covering_read_pairs) %>%
    GRanges()
)




n <- 1

patt <- paste0("chr", n, "[^0-9]")


pytt <- "chr[1-9][0-9]?|chr[XY]"

str_match("xxxxx_chr10.vcf", patt)



lapply(c("chr1", "chr1", "chr5"), function(x){
  GRanges(paste0(x,":1000-5000")) %>%
    return()
}) %>%
  Reduce(function(x,y)c(x,y),.) 







phased <- read_delim("/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-3GNGK9/cnv_results/paired/metastasis02_blood/results_ACEseqWorkflow-1.2.8-4_v1_6_2021-05-22_12h27/phasing/phased_genotype_chr16.txt_haps",
                      col_names = c("id", "start", "end", "ref", "alt", "hp1", "hp2"), delim=" ")


unphased <- read_delim("/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-3GNGK9/cnv_results/paired/metastasis02_blood/results_ACEseqWorkflow-1.2.8-4_v1_6_2021-05-22_12h27/phasing/unphased_chr16.vcf",
                        col_names = c("id", "start", "end", "ref", "alt", "hp1", "hp2"), delim=" ")


hap_block_file <- read_tsv("/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-3GNGK9/cnv_results/paired/metastasis02_blood/results_ACEseqWorkflow-1.2.8-4_v1_6_2021-05-22_12h27/phasing/haploblocks_chr16.txt",
)
### check possible phasing results for POI
found <- c(
"16-23642529" ,   
"16-23642430" ,   
"16-23642600")
lv <- read.vcfR(FILE_FULL_SNV)

View(as_tibble(lv@gt))


comb_v <- bind_cols(as_tibble(lv@fix), as_tibble(lv@gt))


head(comb_v)

chr16 <- comb_v[which(comb_v$CHROM=="16"),]%>%
  
  select(1:9, gtcol=10)


chr16_hq=chr16  %>%
  
  filter(str_detect(gtcol, ":99$"),
         str_detect(gtcol, "255"))


#detailed_genewise_output <- capture.output(

 # , type="message")























#########################
###############################







result <- run_ZP_from_files(files$PID, files$sample, files$seq_method, files$G1_cnv_file, files$rna_file, files$G1_snv_file, files$G1_indel_file,
                            files$G2_file, files$bam_file, raw_vcf, TBL_GENCODE_EXON)






########





## Phasing

#poi <- "H021-2KP9K4_tumor02_blood_WES"
#poi <- "H021-1MA1EB_tumor_buffy_coat_WGS"
poi <- "H021-RGR62Q_tumor_blood_WES"
files <- sf %>% filter(sample_id == poi)


FILE_FULL_SNV <- find_raw_snv(files$PID, files$sample, files$seq_method)
FILE_FULL_INDEL <- find_raw_indel(files$PID, files$sample, files$seq_method)
TBL_GENCODE_EXON <- read_tsv("/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/databases/gencode/gencode19/GencodeV19_Exons_plain.bed.gz",
                             # n,
                             col_names=c("chr", "start", "end","gene", "exon_number", "strand", "width"),
                             skip=1) #%>%
#filter(gene %in% TA_genes)

print("gene_model loaded")

raw_vcf <- c(FILE_FULL_SNV, FILE_FULL_INDEL)

result <- run_ZP_from_files(files$PID, files$sample, files$seq_method, files$G1_cnv_file, files$rna_file, files$G1_snv_file, files$G1_indel_file,
                              files$G2_file, files$bam_file, raw_vcf, TBL_GENCODE_EXON)


out_dir <- file.path("/home/m168r/home_extension/ZP_revision/gt", files$sample_id)
dir.create(out_dir)
store_ZP_results(result, out_dir, files$sample_id)


result$eval_per_variant %>% View()
result$phasing_info %>% View()
result$eval_per_gene








## sub clones


poi <- "H021-9V19Y4"


file <- "/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-9V19Y4/cnv_results/paired/metastasis_blood/results_ACEseqWorkflow-1.2.8-4_v1_6_2023-09-13_06h09/phasing/phased_chr1.vcf"


test <- vcfR::read.vcfR(file)
test@gt %>% as_tibble() %>% View



gr_roi <- GRanges("1:100000-105000")

  tab_vcf <- TabixFile(file)
  vcf_region <- 
    readVcf(tab_vcf, "hg19", 
            param=gr_roi)

#### check if this is also in our germline files
  
  
file <- "/omics/odcf/analysis/hipo/hipo_021/GermlineAnalysis/data_object_master_germline/sequencing/whole_genome_sequencing/results_per_pid/H021-9V19Y4/germline_smallVariants.metastasis.blood/v4.0.5/smallVariants_H021-9V19Y4.clean_annotated.rare.VEP.CharGer.TiNDA.rare_germline.vcf"
  

test <- vcfR::read.vcfR(file)




  
raw_snv_file="/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-D6WTGN/snv_results/paired/metastasis02_blood/results_SNVCallingWorkflow-1.2.166-5_v1_2_2023-06-04_18h52/snvs_H021-D6WTGN.vcf.gz"


test <- vcfR::read.vcfR(raw_snv_file)
  
test@gt %>% as_tibble() %>% .[1:100,]%>% View  
  
  
  






read_tsv(file)

gr_roi <-  GRanges(seqnames = "1", ranges = IRanges(start = 100000, end = 105000))


tab_vcf <- TabixFile(file)


tsv_data <- scanTabix(tab_vcf, which = gr_roi)







vcf_region <- 
  VariantAnnotation::readVcf(tab_vcf, "hg19", 
          param=gr_roi) 

















### ins + del

bamDna <- "/omics/odcf/project/hipo/hipo_021/sequencing/exon_sequencing/view-by-pid/H021-1EEW3X/metastasis/paired/merged-alignment/metastasis_H021-1EEW3X_merged.mdup.bam"
ref_chr1 <- "7"
ref_pos1 <- 4821377
ref_gr1 <- GRanges(seqnames = ref_chr1, 
                   ranges = ref_pos1)
ref_chr2 <- "7"
ref_pos2 <- 4821383
ref_gr2 <- GRanges(seqnames = ref_chr2, 
                   ranges = ref_pos2)
ref_alt1 <- "TG"
ref_alt2 <- "C"
ref_ref1 <- "T"
ref_ref2 <- "CAGGT"
ref_class1 <- "ins"
ref_class2 <- "del"


qname <- "ST-K00246:434:HMJY2BBXY:5:1205:14103:29518"

## now load all reads/read-pairs that cover the position of the first variant
#vm("loading reads", verbose)
all_covering_read_pairs <- readGAlignmentPairs(
  bamDna,
  param=ScanBamParam(
    which=GRanges(seqnames = ref_chr1, 
                  ranges = ref_pos1),
    what=c("qname","seq", "cigar", "mapq", "qual")
  )) 




cigar_string <- "151M"


cigar <- Rsamtools:::.parse_cigar(cigar_string)












all_reads <- c(
  GenomicAlignments::first(all_covering_read_pairs) %>%
    GRanges(),
  GenomicAlignments::last(all_covering_read_pairs) %>%
    GRanges()
)

pileup()

GenomicAlignments::first(all_covering_read_pairs) %>%
  GRanges() %>% as_tibble() %>% View



qual_string <- "--AFJJJAFJJ<JJAFA-JF77-FA<<7-A-AA77---FJA7-FF-A7-"


ascii_encoded <- str_sub(qual_string, start=1, end=1)
exp <- (as.integer(charToRaw(ascii_encoded))-33)/(-10)
# Convert the QUAL string to raw hexadecimal
prob_base_wrong <- 10^exp

# Convert hexadecimal values to numeric Phred quality scores
quality_scores <- as.integer(hex_values) - as.integer(charToRaw('!'))

# Print the quality scores
print(quality_scores)



probs <- c(0.9, 0.95, 0.95,0.95, 0.5)

product(probs)^(1/length(probs))


phased_per_mut <- tibble(mut_id=c("m1", "m2", "m3"), gt=c("1|0", "1|0", "0|1"))

outer(phased_per_mut$mut_id, phased_per_mut$mut_id, `paste`) %>%
  .[which(upper.tri(.))] %>%
  as.data.frame() %>% 
  set_names(nm='raw') %>%
  mutate(
    mut_id1=str_split(raw, " ") %>% map_chr(.,1),
    mut_id2=str_split(raw, " ") %>% map_chr(.,2)
  ) %>%
  left_join(phased_per_mut %>% select(mut_id1=mut_id, gt1=gt), by="mut_id1") %>%
  left_join(phased_per_mut %>% select(mut_id2=mut_id, gt2=gt), by="mut_id2") %>%
  mutate(status=case_when(gt1==gt2 ~ "same",
                          TRUE ~ "diff"),
         comb_id=paste(mut_id1, mut_id2, sep='-')) %>%
  select(comb_id, status)


mut_id1 <- "m3"
mut_id2 <- "m2"

paste(mut_id1, mut_id2, sep='-')
paste(paste0("m",sort(as.numeric(str_match(c(mut_id1, mut_id2), "\\d+")))),
      collapse = "-")






df <- tibble(x=c(1,2,3), y=c(4,5,6), z=c(7,8,9))

add2 <- function(e,r){
  return(e+r)
}

apply(df, 1, add2, e=x, r=y)


for(i in c(1:3)){
  
  print(as.character(x[i,]) %>% set_names(nm=names(x)))
  
}



phasing_dir <- "/omics/odcf/project/hipo/hipo_021/sequencing/whole_genome_sequencing/view-by-pid/H021-3GNGK9/cnv_results/paired/metastasis02_blood/results_ACEseqWorkflow-1.2.8-4_v1_6_2021-05-22_12h27/phasing"


file <- file.path(phasing_dir, "phased_chr1.vcf")
poss <- read.vcfR(file.path(phasing_dir, "phased_chr1.vcf"))
combined <- bind_cols(poss@fix, poss@gt)



pgh <- read_delim(file.path(phasing_dir, "phased_genotype_chr1.txt_haps"),
                  col_names = c("id", "start", "end", "ref", "alt", "hp1", "hp2"), delim=" ")



hap_block_file <-  read_tsv(file.path(phasing_dir, "haploblocks_chr1.txt"), col_names=c("chr", "start", "end", "length"))



gr_roi <- GenomicRanges::GRanges("1:700000-725000")

#indexTabix(file, format = "vcf")

#tab_vcf <- TabixFile(file)
vcf_region <- 
  VariantAnnotation::readVcf(file, "hg19") #%>%
#rowRanges()

all_files <- list.files(phasing_dir, pattern="*")

rem_chr <- str_replace_all(all_files, "chr\\d+", "chrN") %>% str_replace_all("chrX", "chrN") %>%
  unique()


"sample_g.txt"   


"haploblocks_chrN.txt"                     


"phased_chrN.vcf"                          
"phased_genotype_chrN.txt_haps"            
"phased_genotype_chrN.txt_haps_confidence"
"phased_genotype_chrN.txt_info"            
"phased_genotype_chrN.txt_info_by_sample"  
"phased_genotype_chrN.txt_summary"         
"phased_genotype_chrN.txt_warnings"       

"unphased_chrN.vcf"                        
"unphased_genotype_chrN.txt"  













