####载入包####
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
library(ieugwasr)
library(plinkbinr)
library(tidyverse)
library(data.table)
library(rlist)
library(future)
library(MendelianRandomization)
library(TwoSampleMR)
library(MVMR)
library(data.table)
library(MRPRESSO)
library(vcfR)
library(data.table)
library(tidyverse)
library(ggplot2)
rm(list = ls())

####1.暴露数据处理####
setwd("E:/孟德尔自定义函数调试")
source("MR_ALL.R")

#179 lipid#
file_list<-c(paste0("GCST", 90277238:90277268),"GCST90277250","GCST90277288","GCST90277298","GCST90277305","GCST90277311","GCST90277330","GCST90277336")#paste0("GCST", 90277238:90277268),
process_gwas_data(
  work_path="E:/孟德尔自定义函数调试",
  rawdata_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/rawdata", # Path to the directory containing the raw GWAS data files. This is where the function will look for the input files.
  output_no_filter_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/no_filter", # Path to the directory where the function will save the processed data before filtering. This directory will store intermediate results.
  output_clumped_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/clumped", # Path to the directory where the function will save the data after performing clumping. Clumping reduces the number of SNPs by removing those in linkage disequilibrium.
  output_f_filter_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/F_filter", # Path to the directory where the function will save the final filtered data based on the F-statistic. This is the final output of the data processing.
  file_ids = file_list, # A vector of file IDs to process. Each ID corresponds to a GWAS dataset, and the function will process each one in turn.
  file_suffix=".tsv.gz",
  pval_col = "neg_log_10_p_value",
  pFilter=5e-06,
  plink_path = "C:/Users/qifengle/AppData/Local/R/win-library/4.2/plinkbinr/bin/plink_Windows.exe", # Path to the PLINK executable. PLINK is used for clumping SNPs.
  bfile_path = "E:/孟德尔自定义函数调试/gene_ref/EUR/EUR", # Path to the reference genetic data used for clumping. This is typically a binary PLINK file (.bed, .bim, .fam) specifying the reference population.
  metadata_file = "179脂质信息.csv", # Path to the CSV file containing metadata about the GWAS datasets. This file should include sample size and other relevant information.
  metadata_id_col = "accessionId", # The column name in the metadata file that contains the unique IDs for each GWAS dataset. This is used to match the metadata to the GWAS files.
  phenotype_col = "Swiss Lipids Name",
  metadata_sample_size_col = "discoverySampleAncestry" ,# The column name in the metadata file that contains the sample size for each GWAS dataset. This information is used to calculate statistics during processing.
  type= 'outcome', 
  snp_col="rsid",
  beta_col = "beta", 
  se_col = "standard_error", 
  effect_allele_col= "effect_allele" , 
  other_allele_col= "other_allele" , 
  eaf_col= "effect_allele_frequency" , 
  chr_col = "chromosome", 
  pos_col = "base_pair_location"
)



#2.结局数据处理####
rm(list=ls())
#install.packages("vcfR")
#install.packages("geneHapR")
library(vcfR)
library(data.table)
library(tidyverse)
setwd("E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/rawdata")

vcf <- read.vcfR("common_all_20180423.vcf.gz")
vcf2 = vcf@fix[,1:3]
rm(vcf)
fwrite(vcf2,"common_all_20180423.txt",sep = "\t",quote = T,row.names = F,col.names = T)
df<-fread("GCST012879_buildGRCh37.tsv")
colnames(df)
colnames(vcf2)<-c("chromosome","base_pair_location","rsid")
vcf3<-as.data.frame(vcf2)
vcf3$chromosome<-as.numeric(vcf3$chromosome)
vcf3$base_pair_location <-as.numeric(vcf3$base_pair_location)
df_merged <- merge(df, vcf3, by = c('chromosome', 'base_pair_location'), all.x = TRUE)
library(tidyverse)
df_merged2 <- df_merged %>% drop_na("rsid")
df_merged2<-df_merged2[,-c(8:10)]

pFilter_list=c(5e-08)
file<-c("GCST012879")
df<-df_merged2
df$Phenotype<-file
df$n=32072
colnames(df)
df <- format_data(
  df,
  type='outcome',
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col ="effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  phenotype_col = "Phenotype",
  ncase_col = "n",
  chr_col = "chromosome",
  pos_col = "base_pair_location"
)

setwd("E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/no_filter")
fwrite(df,"GCST012879_no_filter.txt",sep = "\t",quote = T,row.names = F,col.names = T)
#读取暴露数据
exposure_dat=read_exposure_data(filename=paste0(file, "_no_filter.txt"),
                                sep= "\t",
                                snp_col = "SNP",
                                beta_col = "beta.outcome",
                                se_col = "se.outcome",
                                pval_col = "pval.outcome",
                                effect_allele_col="effect_allele.outcome",
                                other_allele_col = "other_allele.outcome",
                                eaf_col = "eaf.outcome",
                                phenotype_col = "outcome",
                                samplesize_col = "ncase.outcome",
                                chr_col="chr.outcome", pos_col = "pos.outcome",
                                clump=FALSE)
for (pFilter in pFilter_list) {
  # 筛选数据
  #P值筛选
  exposure_dat_p = subset(exposure_dat, pval.exposure<pFilter)
  # 本地clump
  exposure_dat_clumped <-  ieugwasr::ld_clump(
    # clump_kb = 10000,
    # clump_r2 = 0.001,
    # clump_p = 1,
    dplyr::tibble(rsid=exposure_dat_p $SNP, pval=exposure_dat_p $pval.exposure, id=exposure_dat_p $id.exposure),
    #get_plink_exe()
    plink_bin = "C:/Users/qifengle/AppData/Local/R/win-library/4.2/plinkbinr/bin/plink_Windows.exe",
    #欧洲人群参考基因组位置
    bfile = "E:/孟德尔自定义函数调试/gene_ref/EUR/EUR"
  )
  
  #提取clump后数据
  exposure_dat_clumped<- merge(exposure_dat_p,exposure_dat_clumped, by.x = "SNP",by.y="rsid")

  setwd("E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/clumped/")
  write.table(exposure_dat_clumped, paste0(file,"_",pFilter, "_local_clumped.txt"), sep="\t", row.names=FALSE)
  #计算F值筛选
  exposure_dat_F <- exposure_dat_clumped %>% dplyr::mutate(maf=ifelse(eaf.exposure<0.5,eaf.exposure,(1-eaf.exposure)),
                                                           R2=2*(1-maf)*maf*beta.exposure*beta.exposure,
                                                           F.value=(R2/(1-R2))*((samplesize.exposure-1-1)/1)) %>% dplyr::filter(F.value>10,maf>0.01)
  setwd("E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/F_filter/")
  write.table(exposure_dat_F, paste0(file,"_",pFilter, "_local_clumped_F.txt"), sep="\t", row.names=FALSE)
}






#3.First MR####
perform_mr_analysis(
  exposure_data_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/F_filter", # The directory path where the exposure data files are located. These files contain the genetic variants associated with the exposures.
  outcome_data_path = "E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/no_filter/", # The directory path where the outcome data files are stored. These files contain the genetic associations with the outcome (e.g., colorectal cancer).
  result_output_path = "E:/孟德尔自定义函数调试/MR/", # The directory path where the results of the Mendelian Randomization analysis will be saved. This is where all output files will be written.
  file_list = file_list, # A list of file IDs to process. Each ID corresponds to a specific exposure dataset that will be analyzed in the MR analysis.
  pFilter=5e-06,
  outcome_name = c("GCST012879"), # A vector of outcome names. Each name corresponds to a specific outcome dataset that will be analyzed against the exposures.
  exposure_file_suffix = "_local_clumped_F.txt", # The suffix used for the exposure data files. This suffix helps identify the correct files to be read when processing the exposures.
  outcome_file_suffix="_no_filter.txt"
)



##4.re_mr#######
perform_re_mr_analysis(
  exposure_data_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/F_filter",
  outcome_data_path = "E:/孟德尔自定义函数调试/GWAS_Data/CRC/GCST012879_CRC/no_filter/",
  mr_results_path="E:/孟德尔自定义函数调试/MR/GCST012879",
  result_output_path = "E:/孟德尔自定义函数调试/MR/Re_MR/",
  outcome_name = c("GCST012879"),
  pFilter=5e-06,
  exposure_file_suffix = "_local_clumped_F.txt",
  outcome_file_suffix="_no_filter.txt"
)

