process_gwas_data <- function(rawdata_path, 
                              output_no_filter_path, 
                              output_clumped_path, 
                              output_f_filter_path,
                              work_path,
                              file_ids, 
                              file_suffix , 
                              pFilter, 
                              plink_path, 
                              bfile_path, 
                              metadata_file, 
                              metadata_id_col, 
                              metadata_sample_size_col, 
                              type, # = 'outcome'
                              snp_col, # = "variant_id"
                              beta_col, # = "beta"
                              se_col , #= "standard_error"
                              effect_allele_col , #= "effect_allele"
                              other_allele_col , #= "other_allele"
                              eaf_col , #= "effect_allele_frequency"
                              pval_col, 
                              phenotype_col, 
                              chr_col , #= "chromosome"
                              pos_col ##= "base_pair_location"
                              ) {
  # Load necessary libraries
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
  library(MVMR)
  library(MRPRESSO)
  library(data.table)
  
  # Check if input and output directories exist, create if they don't
  if (!dir.exists(rawdata_path)) {
    stop("Raw data path does not exist.")
  }
  if (!dir.exists(output_no_filter_path)) {
    dir.create(output_no_filter_path, recursive = TRUE)
  }
  if (!dir.exists(output_clumped_path)) {
    dir.create(output_clumped_path, recursive = TRUE)
  }
  if (!dir.exists(output_f_filter_path)) {
    dir.create(output_f_filter_path, recursive = TRUE)
  }
  
  # Read metadata for sample size
  setwd(work_path)
  metadata <- read.csv(metadata_file, header = TRUE, sep = ",", check.names = FALSE)
  
  for (file in file_ids) {
    setwd(rawdata_path)
    filename <- paste0(file, file_suffix)  # Use custom file suffix
    df <- fread(filename)
    df <- df[!is.na(df[[snp_col]]), ]
    
    # Extract SampleSize from metadata using specified columns
    df$SampleSize <- metadata[[metadata_sample_size_col]][metadata[[metadata_id_col]] == file]
    #df$Phenotype <- metadata[[ phenotype_col]][metadata[[metadata_id_col]] == file]
    df$Phenotype<-file
    df$p_value <- 10^(-df[[pval_col]])#
    # Format data
    df <- format_data(
      df,
      type = type,
      snp_col = snp_col,
      beta_col = beta_col,
      se_col = se_col,
      effect_allele_col = effect_allele_col,
      other_allele_col = other_allele_col,
      eaf_col = eaf_col,
      pval_col = "p_value",
      phenotype_col = "Phenotype",
      ncase_col = "SampleSize",
      chr_col = chr_col,
      pos_col = pos_col
    )

    setwd(output_no_filter_path)
    fwrite(df, paste0(file, "_no_filter.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Read exposure data
    setwd(output_no_filter_path)
    exposure_dat <- read_exposure_data(
      filename = paste0(file, "_no_filter.txt"),
      sep = "\t",
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      pval_col = "pval.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      eaf_col = "eaf.outcome",
      phenotype_col = "outcome",
      samplesize_col = "ncase.outcome",
      chr_col = "chr.outcome",
      pos_col = "pos.outcome",
      clump = FALSE
    )
    
    # Filter data by p-value
    exposure_dat_p <- subset(exposure_dat, pval.exposure < pFilter)
    if (nrow(exposure_dat_p) <= 0) {
      next
    }
    
    # Local clumping
    flag <- TRUE
    tryCatch({
      exposure_dat_clumped <- ieugwasr::ld_clump(
        dplyr::tibble(rsid = exposure_dat_p$SNP, pval = exposure_dat_p$pval.exposure, id = exposure_dat_p$id.exposure),
        plink_bin = plink_path,
        bfile = bfile_path
      )
    }, error = function(e) {
      flag <<- FALSE
    })
    if (!flag) next
    
    # Extract data after clumping
    exposure_dat_clumped <- merge(exposure_dat_p, exposure_dat_clumped, by.x = "SNP", by.y = "rsid")
    setwd(output_clumped_path)
    fwrite(exposure_dat_clumped, paste0(file, "_", pFilter, "_local_clumped_F.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Calculate F value for filtering
    exposure_dat_F <- exposure_dat_clumped %>%
      dplyr::mutate(
        maf = ifelse(eaf.exposure < 0.5, eaf.exposure, (1 - eaf.exposure)),
        R2 = 2 * (1 - maf) * maf * beta.exposure * beta.exposure,
        F.value = (R2 / (1 - R2)) * ((samplesize.exposure - 1 - 1) / 1)
      ) %>%
      dplyr::filter(F.value > 10, maf > 0.01)
    
    setwd(output_f_filter_path)
    fwrite(exposure_dat_F, paste0(file, "_", pFilter, "_local_clumped_F.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    print(paste0(file, " is finished!"))
  }
}

# # Example usage:

# process_gwas_data(
#   work_path="E:/孟德尔自定义函数调试",
#   rawdata_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/rawdata", # Path to the directory containing the raw GWAS data files. This is where the function will look for the input files.
#   output_no_filter_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/no_filter", # Path to the directory where the function will save the processed data before filtering. This directory will store intermediate results.
#   output_clumped_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/clumped", # Path to the directory where the function will save the data after performing clumping. Clumping reduces the number of SNPs by removing those in linkage disequilibrium.
#   output_f_filter_path = "E:/孟德尔自定义函数调试/GWAS_Data/179lipid/F_filter", # Path to the directory where the function will save the final filtered data based on the F-statistic. This is the final output of the data processing.
#   file_ids = paste0("GCST", 90277238:90277268), # A vector of file IDs to process. Each ID corresponds to a GWAS dataset, and the function will process each one in turn.
#   file_suffix=".tsv.gz",
#   pval_col = "neg_log_10_p_value",
#   pFilter=5e-06,
#   plink_path = "C:/Users/qifengle/AppData/Local/R/win-library/4.2/plinkbinr/bin/plink_Windows.exe", # Path to the PLINK executable. PLINK is used for clumping SNPs.
#   bfile_path = "E:/孟德尔自定义函数调试/gene_ref/EUR/EUR", # Path to the reference genetic data used for clumping. This is typically a binary PLINK file (.bed, .bim, .fam) specifying the reference population.
#   metadata_file = "179脂质信息.csv", # Path to the CSV file containing metadata about the GWAS datasets. This file should include sample size and other relevant information.
#   metadata_id_col = "accessionId", # The column name in the metadata file that contains the unique IDs for each GWAS dataset. This is used to match the metadata to the GWAS files.
#   phenotype_col = "Swiss Lipids Name",
#   metadata_sample_size_col = "discoverySampleAncestry" ,# The column name in the metadata file that contains the sample size for each GWAS dataset. This information is used to calculate statistics during processing.
#   type= 'outcome', 
#   snp_col="rsid",
#   beta_col = "beta", 
#   se_col = "standard_error", 
#   effect_allele_col= "effect_allele" , 
#   other_allele_col= "other_allele" , 
#   eaf_col= "effect_allele_frequency" , 
#   chr_col = "chromosome", 
#   pos_col = "base_pair_location"
# )


perform_mr_analysis <- function(
    exposure_data_path, 
    outcome_data_path, 
    result_output_path,
    file_list, 
    outcome_name, 
    pFilter,
    exposure_file_suffix, # Customizable file suffix
    outcome_file_suffix
) {
  # Load necessary libraries
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
  
  # Initialize data frames for storing results
  IVW <- data.frame()
  MR_RESULT <- data.frame()
  heter <- data.frame()
  pleio <- data.frame()
  MR_RESULT_IVW <- data.frame()
  presso_Global_result <- data.frame()
  presso_Outlier_result <- data.frame()
  presso_distortion_result <- data.frame()
  presso_result <- data.frame()
  presso_Global_result_cleaned <- data.frame()
  bys <- data.frame()
  cleaned_dat <- data.frame()
  
  
  
  # Check and ensure input directories exist, stop execution if they don't
  if (!dir.exists(exposure_data_path)) {
    stop("Error: Exposure data path does not exist. Please check the path: ", exposure_data_path)
  }
  
  if (!dir.exists(outcome_data_path)) {
    stop("Error: Outcome data path does not exist. Please check the path: ", outcome_data_path)
  }
  
  # Ensure the result output directory exists, create if it doesn't
  if (!dir.exists(result_output_path)) {
    dir.create(result_output_path, recursive = TRUE)
  }
  
  ### 1400 Metabolites_CRC MR Analysis ###
  for (out in outcome_name) {
    outcome_dir <- file.path(result_output_path, out)
    if (!dir.exists(outcome_dir)) {
      dir.create(outcome_dir)
    }
    
    for (file in file_list) {
      setwd(exposure_data_path)
      
      # Generate the file name using the specified pattern
      fileName <- paste0(file, "_", pFilter, exposure_file_suffix)
      
      # Read exposure data
      exposure_dat <- read_exposure_data(
        filename = fileName,
        sep = "\t",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        pval_col = "pval.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure",
        eaf_col = "eaf.exposure",
        phenotype_col = "exposure",
        samplesize_col = "samplesize.exposure",
        chr_col = "chr.exposure", 
        pos_col = "pos.exposure",
        clump = FALSE
      )
      
      outcomeName <- out # Name of the disease displayed in the chart
      
      setwd(outcome_data_path)
      out_file <- paste0(out, outcome_file_suffix)
      
      # Read outcome data with error handling
      flag <- TRUE
      tryCatch({
        outcome_data <- read_outcome_data(
          snps = exposure_dat$SNP,
          filename = out_file, 
          sep = "\t",
          snp_col = "SNP",
          beta_col = "beta.outcome",
          se_col = "se.outcome",
          effect_allele_col = "effect_allele.outcome",
          other_allele_col = "other_allele.outcome",
          pval_col = "pval.outcome",
          eaf_col = "eaf.outcome",
          phenotype_col = "outcome"
        )
      }, error = function(e) {
        flag <<- FALSE
      })
      if (!flag) next
      
      file_dir <- file.path(outcome_dir, file)
      if (!dir.exists(file_dir)) {
        dir.create(file_dir)
      }
      setwd(file_dir)
      
      # Merge exposure and outcome data
      outcome_data$outcome <- outcomeName
      dat <- harmonise_data(exposure_dat, outcome_data)
      
      # Output instruments used for Mendelian randomization
      outTab <- dat[dat$mr_keep == "TRUE", ]
      if (nrow(outTab) <= 3) {
        next
      }
      write.csv(outTab, file = paste0(out, "_", file, "_outcome.csv"), row.names = FALSE)
      
      # MR-PRESSO outlier detection (biased SNP)
      set.seed(56)
      presso <- run_mr_presso(dat)
      write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(file, ".table.MR-PRESSO_Global.csv"))
      write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(file, ".table.MR-PRESSO_Outlier.csv"))
      write.csv(presso[[1]]$`MR-PRESSO results`$`Distortion Test`, file = paste0(file, ".table.MR-PRESSO-distortion-test.csv"))
      
      # Global_test
      mr_presso_Global <- read.csv(paste0(file, ".table.MR-PRESSO_Global.csv"), header = TRUE, sep = ",", check.names = FALSE)
      mr_presso_Global$exposure <- file
      mr_presso_Global$outcome <- out
      presso_Global_result <- rbind(presso_Global_result, mr_presso_Global)
      
      if (mr_presso_Global$Pvalue < 0.05) {
        # Outlier
        mr_presso_Outlier <- read.csv(paste0(file, ".table.MR-PRESSO_Outlier.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_Outlier$exposure <- file
        mr_presso_Outlier$outcome <- out
        mr_presso_Outlier$SNP <- outTab$SNP
        presso_Outlier_result <- rbind(presso_Outlier_result, mr_presso_Outlier)
        
        # Distortion test
        mr_presso_distortion <- read.csv(paste0(file, ".table.MR-PRESSO-distortion-test.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_distortion$exposure <- file
        mr_presso_distortion$outcome <- out
        presso_distortion_result <- rbind(presso_distortion_result, mr_presso_distortion)
        
        outliers_indices <- which(mr_presso_Outlier$Pvalue < 0.05, arr.ind = TRUE)
        outlier_SNPs <- mr_presso_Outlier$SNP[outliers_indices]
        if (length(outlier_SNPs) > 0) {
          dat <- dat[!dat$SNP %in% outlier_SNPs, ]
          write.csv(dat[dat$mr_keep == "TRUE", ], file = paste0(out, "_", file, "_outcome.Cleaned_SNP.csv"), row.names = FALSE)
          outlier_SNP <- data.frame(SNP = outlier_SNPs, exposure = file, outcome = out, Pvaule_Outlier = mr_presso_Outlier$Pvalue[outliers_indices], Pvaule_Global = mr_presso_Global$Pvalue, Pvaule_distortion = mr_presso_distortion$Pvalue)
          presso_result <- rbind(presso_result, outlier_SNP)
        }
      }
      
      cleaned_dat <- rbind(cleaned_dat, dat)
      
      # MR-PRESSO outlier detection (biased SNP)
      if (nrow(dat) > 3) {
        set.seed(56)
        presso <- run_mr_presso(dat)
        write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(file, ".table.MR-PRESSO_Global_cleaned.csv"))
        write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(file, ".table.MR-PRESSO_Outlier_cleaned.csv"))
        write.csv(presso[[1]]$`MR-PRESSO results`$`Distortion Test`, file = paste0(file, ".table.MR-PRESSO-distortion-test_cleaned.csv"))
        
        # Global_test
        mr_presso_Global_cleaned <- read.csv(paste0(file, ".table.MR-PRESSO_Global_cleaned.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_Global_cleaned$exposure <- file
        mr_presso_Global_cleaned$outcome <- out
        presso_Global_result_cleaned <- rbind(presso_Global_result_cleaned, mr_presso_Global_cleaned)
      }
      
      # Mendelian randomization analysis
      mrResult <- mr(dat)
      # Calculate OR values for results
      mrTab <- generate_odds_ratios(mrResult)
      write.csv(mrTab, file = paste0(out, "_VS_", file, "_", "_MRresult.csv"), row.names = FALSE)
      MR_RESULT <- rbind(MR_RESULT, mrTab)
      
      # Mendelian randomization analysis—IVW model
      mrResult_IVW <- mr(dat, method_list = c('mr_ivw_mre', "mr_ivw", "mr_ivw_fe"))
      # Calculate OR values for results—IVW model
      mrTab_IVW <- generate_odds_ratios(mrResult_IVW)
      write.csv(mrTab_IVW, file = paste0(out, "_VS_", file, "_", "IVW_MRresult.csv"), row.names = FALSE)
      MR_RESULT_IVW <- rbind(MR_RESULT_IVW, mrTab_IVW)
      
      # Bayesian weighted model #
      myBWMR <- BWMR(gammahat = dat$beta.exposure,
                     Gammahat = dat$beta.outcome,
                     sigmaX = dat$se.exposure,
                     sigmaY = dat$se.outcome) # run BWMR
      beta <- myBWMR[["beta"]]
      lci95 <- myBWMR[["beta"]] - 1.96 * myBWMR[["se_beta"]]
      uci95 <- myBWMR[["beta"]] + 1.96 * myBWMR[["se_beta"]]
      or <- exp(myBWMR[["beta"]])
      or_lci95 <- exp(lci95)
      or_uci95 <- exp(uci95)
      pval <- myBWMR[["P_value"]]
      myres <- data.frame(method = "BWMR", beta, lci95, uci95, or, or_lci95, or_uci95, pval)
      myres$exposure <- file
      myres$outcome <- out
      write.csv(myres, file = paste0(out, "_VS_", file, "_", "BWMR.csv"), row.names = FALSE)
      bys <- rbind(bys, myres)
      
      # Heterogeneity analysis
      heterTab <- mr_heterogeneity(dat)
      write.csv(heterTab, file = paste0(out, "_VS_", file, "_", "_heterogeneity.csv"), row.names = FALSE)
      heter <- rbind(heter, heterTab)
      
      # Pleiotropy test
      pleioTab <- mr_pleiotropy_test(dat)
      write.csv(pleioTab, file = paste0(out, "_VS_", file, "_", "_pleiotropy.csv"), row.names = FALSE)
      pleio <- rbind(pleio, pleioTab)
      
      # Generate scatter plot
      pdf(file = paste0(file, "_VS_", out, "_", ".scatter_plot.pdf"), width = 7, height = 6.5)
      p1 <- mr_scatter_plot(mrResult, dat)
      print(p1)
      dev.off()
      
      # Forest plot
      res_single <- mr_singlesnp(dat)      # Obtain the effect of each instrument variable on the outcome
      pdf(file = paste0(file, "_VS_", out, "_", ".forest.pdf"), width = 6.5, height = 5)
      p2 <- mr_forest_plot(res_single)
      print(p2)
      dev.off()
      
      # Funnel plot
      pdf(file = paste0(file, "_VS_", out, "_", ".funnel_plot.pdf"), width = 6.5, height = 6)
      p3 <- mr_funnel_plot(singlesnp_results = res_single)
      print(p3)
      dev.off()
      
      # Leave-one-out sensitivity analysis
      pdf(file = paste0(file, "_VS_", out, "_", ".leaveoneout.pdf"), width = 6.5, height = 5)
      p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
      print(p4)
      dev.off()
      
      #### Summarize IVW results
      mrFile <- paste0(out, "_VS_", file, "_", "_MRresult.csv")      # Mendelian randomization analysis results file
      pleFile <- paste0(out, "_VS_", file, "_", "_pleiotropy.csv")     # Pleiotropy results file
      
      # Read Mendelian randomization results file
      rt <- read.csv(mrFile, header = TRUE, sep = ",", check.names = FALSE)
      # Extract IVW method pvalue < 0.05
      ivw <- rt[((rt$method == "Inverse variance weighted") & (rt$pval < 0.05)), ]
      
      # Read pleiotropy results file
      pleRT <- read.csv(pleFile, header = TRUE, sep = ",", check.names = FALSE)
      # Exclude pvalue less than 0.05, pleiotropy not allowed
      pleRT <- pleRT[pleRT$pval > 0.05, ]
      gutLists <- as.vector(pleRT$exposure)
      outTab <- ivw[ivw$exposure %in% gutLists, ]
      IVW <- rbind(IVW, outTab)
      write.csv(outTab, file = paste0(out, "_VS_", file, "_", "IVW.filter.csv"), row.names = FALSE)
      
      # If the number of rows in exposure_dat is less than or equal to 3, skip the current loop and continue with the next loop
      if (nrow(outTab) == 0) {
        next
      }
      print(paste0(file, " is finished!"))
    }
    
    setwd(outcome_dir)
    IVW_FDR <- subset(MR_RESULT, method == 'Inverse variance weighted')
    IVW_FDR$P.adj <- p.adjust(IVW_FDR$pval, method = "fdr")
    write.csv(IVW_FDR, file = paste0(out, "_MR_Result_FDR_ALL.CSV"))
    
    save.image(file = paste0(out, "_Result_ALL.RData"))
    write.csv(MR_RESULT, file = paste0(out, "_MR_Result_ALL.CSV"))
    write.csv(heter, file = paste0(out, "_MR_heter_ALL.CSV"))
    write.csv(pleio, file = paste0(out, "_MR_pleio_ALL.CSV"))
    write.csv(IVW, file = paste0(out, "_MR_IVW_ALL.CSV"))
    write.csv(MR_RESULT_IVW, file = paste0(out, "_MR_IVW_model_ALL.CSV"))
    write.csv(presso_Outlier_result, file = paste0(out, "_MR_presso_Outlier_result_ALL.CSV"))
    write.csv(presso_distortion_result, file = paste0(out, "_MR_presso_distortion_result_ALL.CSV"))
    write.csv(presso_Global_result, file = paste0(out, "_MR_presso_Global_result_ALL.CSV"))
    write.csv(presso_result, file = paste0(out, "_MR_presso_result_ALL.CSV"))
    write.csv(presso_Global_result_cleaned, file = paste0(out, "_MR_presso_Global_result_cleaned_ALL.CSV"))
    write.csv(bys, file = paste0(out, "_BWMR_ALL.CSV"))
    write.csv(cleaned_dat, file = paste0(out, "_cleaned_outcome_snp_ALL.CSV"))
    
    # Clear data frames for next outcome
    rm(MR_RESULT)
    rm(heter)
    rm(pleio)
    rm(IVW)
    rm(MR_RESULT_IVW)
    rm(presso_Outlier_result)
    rm(presso_distortion_result)
    rm(presso_Global_result)
    rm(bys)
    rm(presso_result)
    rm(cleaned_dat)
  }
}

# Example usage:
# perform_mr_analysis(
#   exposure_data_path = "E:/test/F_filter/", # The directory path where the exposure data files are located. These files contain the genetic variants associated with the exposures.
#   outcome_data_path = "F:/孟德尔随机化--数据+分析/孟德尔数据/结直肠癌/Finngen/5e-08", # The directory path where the outcome data files are stored. These files contain the genetic associations with the outcome (e.g., colorectal cancer).
#   result_output_path = "E:/test/MR/", # The directory path where the results of the Mendelian Randomization analysis will be saved. This is where all output files will be written.
#   BWMR_PATH = "E:/test/", # The directory path where additional necessary scripts (like BWMR_updated.R) are located. This script is sourced for additional computations required during the MR analysis.
#   file_list = paste0("GCST90", 199621:199624), # A list of file IDs to process. Each ID corresponds to a specific exposure dataset that will be analyzed in the MR analysis.
#   outcome_name = c("Finngen_R10_colorectal_cancer"), # A vector of outcome names. Each name corresponds to a specific outcome dataset that will be analyzed against the exposures.
#   exposure_file_suffix = "_local_clumped_F.txt" # The suffix used for the exposure data files. This suffix helps identify the correct files to be read when processing the exposures.
# )


perform_re_mr_analysis <- function(
    exposure_data_path, 
    outcome_data_path, 
    result_output_path,
    mr_results_path,
    outcome_name, 
    pFilter,
    exposure_file_suffix, # Customizable file suffix
    outcome_file_suffix
) {
  # Load necessary libraries
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
  
  # Initialize data frames for storing results
  IVW <- data.frame()
  MR_RESULT <- data.frame()
  heter <- data.frame()
  pleio <- data.frame()
  MR_RESULT_IVW <- data.frame()
  presso_Global_result <- data.frame()
  presso_Outlier_result <- data.frame()
  presso_distortion_result <- data.frame()
  presso_result <- data.frame()
  presso_Global_result_cleaned <- data.frame()
  bys <- data.frame()
  cleaned_dat <- data.frame()
  
  # Check and ensure input directories exist, stop execution if they don't
  if (!dir.exists(exposure_data_path)) {
    stop("Error: Exposure data path does not exist. Please check the path: ", exposure_data_path)
  }
  
  if (!dir.exists(outcome_data_path)) {
    stop("Error: Outcome data path does not exist. Please check the path: ", outcome_data_path)
  }
  if (!dir.exists(mr_results_path)) {
    stop("Error: Outcome data path does not exist. Please check the path: ", mr_results_path)
  }
  # Ensure the result output directory exists, create if it doesn't
  if (!dir.exists(result_output_path)) {
    dir.create(result_output_path, recursive = TRUE)
  }
  
  
  for (out in outcome_name) {
    outcome_dir <- file.path(result_output_path, out)
    if (!dir.exists(outcome_dir)) {
      dir.create(outcome_dir)
    }
    
    # Retrieve file_list from IVW results
    setwd(mr_results_path)
    ivw_results_file <- paste0(out, "_MR_IVW_ALL.CSV")
    rt <- read.csv(ivw_results_file, header = TRUE, sep = ",", check.names = FALSE)
    file_list <- rt$exposure
    
    for (file in file_list) {
      setwd(exposure_data_path)
      
      # Generate the file name using the specified pattern
      fileName <- paste0(file, "_", pFilter, exposure_file_suffix)
      
      # Read exposure data
      exposure_dat <- read_exposure_data(
        filename = fileName,
        sep = "\t",
        snp_col = "SNP",
        beta_col = "beta.exposure",
        se_col = "se.exposure",
        pval_col = "pval.exposure",
        effect_allele_col = "effect_allele.exposure",
        other_allele_col = "other_allele.exposure",
        eaf_col = "eaf.exposure",
        phenotype_col = "exposure",
        samplesize_col = "samplesize.exposure",
        chr_col = "chr.exposure", 
        pos_col = "pos.exposure",
        clump = FALSE
      )
      
      outcomeName <- out # Name of the disease displayed in the chart
      
      setwd(outcome_data_path)
      out_file <- paste0(out, outcome_file_suffix)
      
      # Read outcome data with error handling
      flag <- TRUE
      tryCatch({
        outcome_data <- read_outcome_data(
          snps = exposure_dat$SNP,
          filename = out_file, 
          sep = "\t",
          snp_col = "SNP",
          beta_col = "beta.outcome",
          se_col = "se.outcome",
          effect_allele_col = "effect_allele.outcome",
          other_allele_col = "other_allele.outcome",
          pval_col = "pval.outcome",
          eaf_col = "eaf.outcome",
          phenotype_col = "outcome"
        )
      }, error = function(e) {
        message("Error reading outcome data for ", file, ": ", e$message)
        flag <<- FALSE
      })
      if (!flag) next
      
      file_dir <- file.path(outcome_dir, file)
      if (!dir.exists(file_dir)) {
        dir.create(file_dir)
      }
      setwd(file_dir)
      
      # Merge exposure and outcome data
      outcome_data$outcome <- outcomeName
      dat <- harmonise_data(exposure_dat, outcome_data)
      
      # Output instruments used for Mendelian randomization
      outTab <- dat[dat$mr_keep == "TRUE", ]
      if (nrow(outTab) <= 3) {
        next
      }
      write.csv(outTab, file = paste0(out, "_", file, "_outcome.csv"), row.names = FALSE)
      
      # MR-PRESSO outlier detection (biased SNP)
      set.seed(56)
      presso <- run_mr_presso(dat)
      write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(file, ".table.MR-PRESSO_Global.csv"))
      write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(file, ".table.MR-PRESSO_Outlier.csv"))
      write.csv(presso[[1]]$`MR-PRESSO results`$`Distortion Test`, file = paste0(file, ".table.MR-PRESSO-distortion-test.csv"))
      
      # Global_test
      mr_presso_Global <- read.csv(paste0(file, ".table.MR-PRESSO_Global.csv"), header = TRUE, sep = ",", check.names = FALSE)
      mr_presso_Global$exposure <- file
      mr_presso_Global$outcome <- out
      presso_Global_result <- rbind(presso_Global_result, mr_presso_Global)
      
      if (mr_presso_Global$Pvalue < 0.05) {
        # Outlier
        mr_presso_Outlier <- read.csv(paste0(file, ".table.MR-PRESSO_Outlier.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_Outlier$exposure <- file
        mr_presso_Outlier$outcome <- out
        mr_presso_Outlier$SNP <- outTab$SNP
        presso_Outlier_result <- rbind(presso_Outlier_result, mr_presso_Outlier)
        
        # Distortion test
        mr_presso_distortion <- read.csv(paste0(file, ".table.MR-PRESSO-distortion-test.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_distortion$exposure <- file
        mr_presso_distortion$outcome <- out
        presso_distortion_result <- rbind(presso_distortion_result, mr_presso_distortion)
        
        outliers_indices <- which(mr_presso_Outlier$Pvalue < 0.05, arr.ind = TRUE)
        outlier_SNPs <- mr_presso_Outlier$SNP[outliers_indices]
        if (length(outlier_SNPs) > 0) {
          dat <- dat[!dat$SNP %in% outlier_SNPs, ]
          write.csv(dat[dat$mr_keep == "TRUE", ], file = paste0(out, "_", file, "_outcome.Cleaned_SNP.csv"), row.names = FALSE)
          outlier_SNP <- data.frame(SNP = outlier_SNPs, exposure = file, outcome = out, Pvaule_Outlier = mr_presso_Outlier$Pvalue[outliers_indices], Pvaule_Global = mr_presso_Global$Pvalue, Pvaule_distortion = mr_presso_distortion$Pvalue)
          presso_result <- rbind(presso_result, outlier_SNP)
        }
      }
      
      cleaned_dat <- rbind(cleaned_dat, dat)
      
      # MR-PRESSO outlier detection (biased SNP)
      if (nrow(dat) > 3) {
        set.seed(56)
        presso <- run_mr_presso(dat)
        write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = paste0(file, ".table.MR-PRESSO_Global_cleaned.csv"))
        write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = paste0(file, ".table.MR-PRESSO_Outlier_cleaned.csv"))
        write.csv(presso[[1]]$`MR-PRESSO results`$`Distortion Test`, file = paste0(file, ".table.MR-PRESSO-distortion-test_cleaned.csv"))
        
        # Global_test
        mr_presso_Global_cleaned <- read.csv(paste0(file, ".table.MR-PRESSO_Global_cleaned.csv"), header = TRUE, sep = ",", check.names = FALSE)
        mr_presso_Global_cleaned$exposure <- file
        mr_presso_Global_cleaned$outcome <- out
        presso_Global_result_cleaned <- rbind(presso_Global_result_cleaned, mr_presso_Global_cleaned)
      }
      
      # Mendelian randomization analysis
      mrResult <- mr(dat)
      # Calculate OR values for results
      mrTab <- generate_odds_ratios(mrResult)
      write.csv(mrTab, file = paste0(out, "_VS_", file, "_", "_MRresult.csv"), row.names = FALSE)
      MR_RESULT <- rbind(MR_RESULT, mrTab)
      
      # Mendelian randomization analysis—IVW model
      mrResult_IVW <- mr(dat, method_list = c('mr_ivw_mre', "mr_ivw", "mr_ivw_fe"))
      # Calculate OR values for results—IVW model
      mrTab_IVW <- generate_odds_ratios(mrResult_IVW)
      write.csv(mrTab_IVW, file = paste0(out, "_VS_", file, "_", "IVW_MRresult.csv"), row.names = FALSE)
      MR_RESULT_IVW <- rbind(MR_RESULT_IVW, mrTab_IVW)
      
      # Bayesian weighted model #
      myBWMR <- BWMR(gammahat = dat$beta.exposure,
                     Gammahat = dat$beta.outcome,
                     sigmaX = dat$se.exposure,
                     sigmaY = dat$se.outcome) # run BWMR
      beta <- myBWMR[["beta"]]
      lci95 <- myBWMR[["beta"]] - 1.96 * myBWMR[["se_beta"]]
      uci95 <- myBWMR[["beta"]] + 1.96 * myBWMR[["se_beta"]]
      or <- exp(myBWMR[["beta"]])
      or_lci95 <- exp(lci95)
      or_uci95 <- exp(uci95)
      pval <- myBWMR[["P_value"]]
      myres <- data.frame(method = "BWMR", beta, lci95, uci95, or, or_lci95, or_uci95, pval)
      myres$exposure <- file
      myres$outcome <- out
      write.csv(myres, file = paste0(out, "_VS_", file, "_", "BWMR.csv"), row.names = FALSE)
      bys <- rbind(bys, myres)
      
      # Heterogeneity analysis
      heterTab <- mr_heterogeneity(dat)
      write.csv(heterTab, file = paste0(out, "_VS_", file, "_", "_heterogeneity.csv"), row.names = FALSE)
      heter <- rbind(heter, heterTab)
      
      # Pleiotropy test
      pleioTab <- mr_pleiotropy_test(dat)
      write.csv(pleioTab, file = paste0(out, "_VS_", file, "_", "_pleiotropy.csv"), row.names = FALSE)
      pleio <- rbind(pleio, pleioTab)
      
      # Generate scatter plot
      pdf(file = paste0(file, "_VS_", out, "_", ".scatter_plot.pdf"), width = 7, height = 6.5)
      p1 <- mr_scatter_plot(mrResult, dat)
      print(p1)
      dev.off()
      
      # Forest plot
      res_single <- mr_singlesnp(dat)      # Obtain the effect of each instrument variable on the outcome
      pdf(file = paste0(file, "_VS_", out, "_", ".forest.pdf"), width = 6.5, height = 5)
      p2 <- mr_forest_plot(res_single)
      print(p2)
      dev.off()
      
      # Funnel plot
      pdf(file = paste0(file, "_VS_", out, "_", ".funnel_plot.pdf"), width = 6.5, height = 6)
      p3 <- mr_funnel_plot(singlesnp_results = res_single)
      print(p3)
      dev.off()
      
      # Leave-one-out sensitivity analysis
      pdf(file = paste0(file, "_VS_", out, "_", ".leaveoneout.pdf"), width = 6.5, height = 5)
      p4 <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
      print(p4)
      dev.off()
      
      #### Summarize IVW results
      mrFile <- paste0(out, "_VS_", file, "_", "_MRresult.csv")      # Mendelian randomization analysis results file
      pleFile <- paste0(out, "_VS_", file, "_", "_pleiotropy.csv")     # Pleiotropy results file
      
      # Read Mendelian randomization results file
      rt <- read.csv(mrFile, header = TRUE, sep = ",", check.names = FALSE)
      # Extract IVW method pvalue < 0.05
      ivw <- rt[((rt$method == "Inverse variance weighted") & (rt$pval < 0.05)), ]
      
      # Read pleiotropy results file
      pleRT <- read.csv(pleFile, header = TRUE, sep = ",", check.names = FALSE)
      # Exclude pvalue less than 0.05, pleiotropy not allowed
      pleRT <- pleRT[pleRT$pval > 0.05, ]
      gutLists <- as.vector(pleRT$exposure)
      outTab <- ivw[ivw$exposure %in% gutLists, ]
      IVW <- rbind(IVW, outTab)
      write.csv(outTab, file = paste0(out, "_VS_", file, "_", "IVW.filter.csv"), row.names = FALSE)
      
      # If the number of rows in exposure_dat is less than or equal to 3, skip the current loop and continue with the next loop
      if (nrow(outTab) == 0) {
        next
      }
      print(paste0(file, " is finished!"))
    }
    
    setwd(outcome_dir)
    save.image(file = paste0(out, "_Result.RData"))
    write.csv(MR_RESULT, file = paste0(out, "_MR_Result.CSV"))
    write.csv(heter, file = paste0(out, "_MR_heter.CSV"))
    write.csv(pleio, file = paste0(out, "_MR_pleio.CSV"))
    write.csv(IVW, file = paste0(out, "_MR_IVW.CSV"))
    write.csv(MR_RESULT_IVW, file = paste0(out, "_MR_IVW_model.CSV"))
    write.csv(presso_Outlier_result, file = paste0(out, "_MR_presso_Outlier_result.CSV"))
    write.csv(presso_distortion_result, file = paste0(out, "_MR_presso_distortion_result.CSV"))
    write.csv(presso_Global_result, file = paste0(out, "_MR_presso_Global_result.CSV"))
    write.csv(presso_result, file = paste0(out, "_MR_presso_result.CSV"))
    write.csv(presso_Global_result_cleaned, file = paste0(out, "_MR_presso_Global_result_cleaned.CSV"))
    write.csv(bys, file = paste0(out, "_BWMR.CSV"))
    write.csv(cleaned_dat, file = paste0(out, "_cleaned_outcome_snp.CSV"))
    
    # Clear data frames for next outcome
    rm(MR_RESULT, heter, pleio, IVW, MR_RESULT_IVW, presso_Outlier_result, presso_distortion_result, presso_Global_result, bys, presso_result, cleaned_dat)
  }
}

# # Example usage:
# perform_re_mr_analysis(
#   exposure_data_path = "F:/孟德尔随机化--数据+分析/孟德尔数据/1400代谢物/F_filter/5e-06/",
#   outcome_data_path = "F:/孟德尔随机化--数据+分析/孟德尔数据/结直肠癌/Finngen/5e-08",
#   mr_results_path="F:/孟德尔随机化--数据+分析/孟德尔分析结果汇总/1400代谢物-免疫表型-肠道菌群/1400 metabolites and colorectal cancer/1400代谢物_CRC_MR_outcome_5e-06/Finngen_R10_colorectal_cancer",
#   result_output_path = "E:/test/MR/Re_MR/",
#   BWMR_PATH = "E:/test/",
#   outcome_name = c("Finngen_R10_colorectal_cancer"),
#   pFilter=5e-06,
#   exposure_file_suffix = "local_clumped_F.txt"
# )

# BWMR_updated.R content begins
###### VEM and Inference for BWMR ######
##
### INPUT
## gammahat:                        SNP-exposure effect;
## Gammahat:                        SNP-outcome effect;
## sigmaX:                          standard error of SNP-exposure effect;
## sigmaY:                          standard error of SNP-outcome effect;
##
### OUTPUT
## beta                             the estimate of beta;
## se_beta                          the estimate the standard error of beta;
## P-value                          P-value
## plot1                            Plot of Data with Standard Error Bar;
## plot2                            Plot of Evidence Lower Bound (ELBO);
## plot3                            Posterior Mean of Weight of Each Observation;
## plot4                            Plot of Weighted Data and Its Regression Result.


# packages
library("ggplot2")


# known parameters for the prior distributions
sqsigma0 <- (1e+6)^2
alpha <- 100


## define function to calculate ELBO and E[Lc] (the approximate log-likelihood)
ELBO_func <- function(N, gammahat, Gammahat, sqsigmaX, sqsigmaY, mu_beta, sqsigma_beta, mu_gamma, sqsigma_gamma, a, b, pi_w, sqsigma, sqtau) {
  # + E[log p(gammahat, sqsigmaX | gamma)]
  l <- - 0.5*sum(log(sqsigmaX)) - 0.5*sum(((gammahat-mu_gamma)^2+sqsigma_gamma)/sqsigmaX)    
  # + E[log p(Gammahat, sqsigmaY| beta, gamma, w, sqtau)]
  l <- l - 0.5*log(2*pi)*sum(pi_w) - 0.5*sum(pi_w*log(sqsigmaY+sqtau)) 
  l <- l - 0.5*sum(pi_w*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)/(sqsigmaY+sqtau))
  # + E[log p(beta | sqsigma0)]
  l <- l - 0.5*(mu_beta^2+sqsigma_beta)/sqsigma0
  # + E[log p(gamma | sqsigma)]
  l <- l - 0.5*N*log(sqsigma) - 0.5/sqsigma*sum(mu_gamma^2+sqsigma_gamma)
  # + E[log p(w | pi1)]
  l <- l + (digamma(a)-digamma(a+b))*sum(pi_w) + (digamma(b)-digamma(a+b))*(N-sum(pi_w))
  # + E[log p(pi1)]
  l <- l + (alpha-1)*(digamma(a)-digamma(a+b))
  
  # - E[log q(beta)]
  l <- l + 0.5*log(sqsigma_beta)
  # - E[log q(gamma)]
  l <- l + 0.5*sum(log(sqsigma_gamma))
  # - E[log q(pi1)]
  l <- l - (a-1)*(digamma(a)-digamma(a+b)) - (b-1)*(digamma(b)-digamma(a+b)) + lbeta(a, b)
  # - E[log q(w)]
  # Need to check if pi_w = 0 or pi_w = 1, since there are log terms of pi_w and 1 - pi_w.
  # e1 <- pi_w*log(pi_w)
  # e2 <- (1-pi_w)*log(1-pi_w)
  # e1[which(pi_w == 0)] <- 0
  # e2[which(pi_w == 1)] <- 0
  # l <- l - sum(e1+e2)
  ## A TRICK TO SIMPLIFY THE CODE
  l <- l - sum(pi_w*log(pi_w+(pi_w==0)) + (1-pi_w)*log(1-pi_w+(pi_w==1)))
  # l: ELBO
}


BWMR <- function(gammahat, Gammahat, sigmaX, sigmaY) {
  ## data
  N <- length(gammahat)
  sqsigmaX <- sigmaX^2
  sqsigmaY <- sigmaY^2
  
  ### Variational EM algorithm ###
  # initialize
  # initial parameters of BWMR
  beta <- 0
  sqtau <- 1^2
  sqsigma <- 1^2
  # initial parameters of variational distribution
  mu_gamma <- gammahat
  sqsigma_gamma <- rep(0.1, N)
  pi_w <- rep(0.5, N)
  # declare sets of ELBO and approximate log-likelihood
  ELBO_set <- numeric(0)
  
  for (iter in 1:5000) {
    ## Variational E-Step
    # beta
    sqsigma_beta <- 1/(1/sqsigma0 + sum(pi_w*(mu_gamma^2+sqsigma_gamma)/(sqsigmaY+sqtau)))
    mu_beta <- sum(pi_w*mu_gamma*Gammahat/(sqsigmaY+sqtau))*sqsigma_beta
    # gamma
    sqsigma_gamma <- 1/(1/sqsigmaX + pi_w*(mu_beta^2+sqsigma_beta)/(sqsigmaY+sqtau) + 1/sqsigma)
    mu_gamma <- (gammahat/sqsigmaX + pi_w*Gammahat*mu_beta/(sqsigmaY+sqtau))*sqsigma_gamma
    # pi1
    a <- alpha + sum(pi_w)
    b <- N + 1 - sum(pi_w)
    # w
    q0 <- exp(digamma(b) - digamma(a+b))
    q1 <- exp(- 0.5*log(2*pi) - 0.5*log(sqsigmaY+sqtau) - 0.5*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)/(sqsigmaY+sqtau) + digamma(a)-digamma(a+b))
    pi_w <- q1/(q0+q1)
    
    if (sum(pi_w) == 0){
      message("Invalid IVs!")
      mu_beta = NA
      se_beta = NA
      P_value = NA
      return(list(beta=NA, se_beta=NA, P_value=NA))
    }
    
    ## M-Step
    sqsigma <- sum(mu_gamma^2 + sqsigma_gamma)/N
    sqtau <- sum(pi_w*((mu_beta^2+sqsigma_beta)*(mu_gamma^2+sqsigma_gamma)-2*mu_beta*mu_gamma*Gammahat+Gammahat^2)*sqtau^2/((sqsigmaY+sqtau)^2)) / sum(pi_w/(sqsigmaY+sqtau))
    sqtau <- sqrt(sqtau)
    
    ## check ELBO
    ELBO <- ELBO_func(N, gammahat, Gammahat, sqsigmaX, sqsigmaY, mu_beta, sqsigma_beta, mu_gamma, sqsigma_gamma, a, b, pi_w, sqsigma, sqtau)
    ELBO_set <- c(ELBO_set, ELBO)
    if (iter > 1 && (abs((ELBO_set[iter]-ELBO_set[iter-1])/ELBO_set[iter-1]) < 1e-6)) {
      break
    }
  }
  # message("Iteration=", iter, ", beta=", mu_beta, ", tau=", sqrt(sqtau), ", sigma=", sqrt(sqsigma), ".")
  
  
  
  ### visualize the result of VEM algorithm
  # Plot1: Plot of Data with Standard Error Bar
  df1 <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sigmaX = sigmaX,
    sigmaY = sigmaY
  )
  plot1 <- ggplot(data = df1, aes(x = gammahat, y = Gammahat)) +  
    geom_pointrange(aes(ymin = Gammahat - sigmaY, ymax = Gammahat + sigmaY), color="gray59", size = 0.3) +
    geom_errorbarh(aes(xmin = gammahat - sigmaX, xmax = gammahat + sigmaX, height = 0), color="gray59") +
    labs(x = "SNP-exposure effect", y = "SNP-outcome effect", title = "Plot1: Plot of data with standard error bar")
  
  # Plot2: Plot of Evidence Lower Bound (ELBO)
  iteration <- seq(1, (length(ELBO_set)), by = 1)
  df2 <- data.frame(
    iteration = iteration,
    ELBO_iter = ELBO_set
  )
  plot2 <- ggplot(df2, aes(x=iteration, y=ELBO_iter)) + geom_line(size = 0.5, color = "tomato1") + geom_point(size=0.5, color = "tomato1") +
    labs(x = "iteration", y="elbo", title = "Plot2: Plot of evidence lower bound (elbo)")
  
  # Plot3: Posterior Mean of Weight of Each Observation
  serial_number <- seq(1, N, by = 1)
  df3 <- data.frame(
    weight = pi_w,
    serial_number = serial_number
  )
  plot3 <- ggplot(data = df3, mapping = aes(x = factor(serial_number), y = weight, fill = weight)) + geom_bar(stat = 'identity', position = 'dodge') +
    labs(x = "observation No.", y = "weight", title = "Plot3: Posterior mean of weight of each observation") +
    ylim(0, 1) +
    theme(axis.text.x = element_text(size = 5))
  # scale_x_discrete(breaks = seq(10, N, 20)) +
  
  # Plot4: Plot of Weighted Data and Its Regression Result
  df4 <- data.frame(
    gammahat = gammahat,
    Gammahat = Gammahat,
    sqsigmaX = sqsigmaX,
    sqsigmaY = sqsigmaY,
    w = pi_w
  )
  plot4 <- ggplot(df4, aes(x=gammahat, y=Gammahat, color=w)) + geom_point(size = 0.3) +
    geom_pointrange(aes(ymin = Gammahat - sigmaY, ymax = Gammahat + sigmaY), size = 0.3) +
    geom_errorbarh(aes(xmin = gammahat - sigmaX, xmax = gammahat + sigmaX, height = 0)) +
    geom_abline(intercept=0, slope=mu_beta, color="#990000", linetype="dashed", size=0.5) +
    labs(x = "SNP-exposure effect", y = "SNP-outcome effect", title = "Plot4: Plot of weighted data and its regression result")
  
  
  ### LRVB and Standard Error ###
  ## matrix V
  forV <- matrix(nrow = N, ncol = 4)
  forV[ ,1] <- sqsigma_gamma
  forV[ ,2] <- 2*mu_gamma*sqsigma_gamma
  forV[ ,3] <- forV[ ,2]
  forV[ ,4] <- 2*sqsigma_gamma^2 + 4*mu_gamma^2*sqsigma_gamma
  V <- matrix(rep(0, (3*N+4)*(3*N+4)), nrow = 3*N+4, ncol = 3*N+4)
  for (j in 1:N) {
    V[(3*j):(3*j+1), (3*j):(3*j+1)] <- matrix(forV[j, ], 2, 2)
    V[3*j+2, 3*j+2] <- pi_w[j] - (pi_w[j]^2)
  }
  V[1:2, 1:2] <- matrix(c(sqsigma_beta, 2*mu_beta*sqsigma_beta, 2*mu_beta*sqsigma_beta, 2*sqsigma_beta^2+4*mu_beta^2*sqsigma_beta), 2, 2)
  V[(3*N+3):(3*N+4), (3*N+3):(3*N+4)] <- matrix(c(trigamma(a)-trigamma(a+b), -trigamma(a+b), -trigamma(a+b), trigamma(b)-trigamma(a+b)), 2, 2)
  
  ## matrix H
  H <- matrix(rep(0, (3*N+4)*(3*N+4)), nrow = 3*N+4, ncol = 3*N+4)
  forH <- matrix(nrow = N, ncol = 6)
  forH[ ,1] <- pi_w*Gammahat/(sqsigmaY+sqtau)
  forH[ ,2] <- mu_gamma*Gammahat/(sqsigmaY+sqtau)
  forH[ ,3] <- -0.5*pi_w/(sqsigmaY+sqtau)
  forH[ ,4] <- -0.5*(mu_gamma^2+sqsigma_gamma)/(sqsigmaY+sqtau)
  forH[ ,5] <- mu_beta*Gammahat/(sqsigmaY+sqtau)
  forH[ ,6] <- -0.5*(mu_beta^2+sqsigma_beta)/(sqsigmaY+sqtau)
  for (j in 1:N) {
    H[1, 3*j] <- forH[j, 1]
    H[1, 3*j+2] <- forH[j, 2]
    H[2, 3*j+1] <- forH[j, 3]
    H[2, 3*j+2] <- forH[j, 4]
    H[(3*N+3):(3*N+4), 3*j+2] <- c(1, -1)
    H[3*j+2, (3*N+3):(3*N+4)] <- c(1, -1)
    H[(3*j):(3*j+1), 3*j+2] <- forH[j, 5:6]
    H[3*j+2, (3*j):(3*j+1)] <- forH[j, 5:6]
  }
  H[ ,1] <- H[1, ]
  H[ ,2] <- H[2, ]
  
  
  ## accurate covariance estimate and standard error
  I <- diag(3*N+4)
  Sigma_hat <- try(solve(I-V%*%H)%*%V)
  
  if (inherits(Sigma_hat, "try-error")){
    message("Invalid IVs!")
    return(list(beta=NA, se_beta=NA, P_value=NA))
  } else{
    se_beta <- sqrt(Sigma_hat[1, 1])
  }
  
  
  ## test
  W <- (mu_beta/se_beta)^2
  # P_value <- 1 - pchisq(W, 1)
  P_value <- pchisq(W, 1, lower.tail=F)
  
  message("Estimate of beta=", mu_beta, ", se of beta=", se_beta, ", P-value=", P_value, ".")
  
  ## output
  output <- list(beta=mu_beta, se_beta=se_beta, P_value=P_value, 
                 weights=pi_w, tau=sqrt(sqtau), sigma=sqrt(sqsigma), mu_pi=a/(a+b),
                 plot1=plot1, plot2=plot2, plot3=plot3, plot4=plot4)
}
