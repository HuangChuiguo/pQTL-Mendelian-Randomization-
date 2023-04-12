#######################################################################################################################################
#Two-sample Mendelian randomization 
#######################################################################################################################################
rm(list = ls())
library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(tidyverse)
library(TwoSampleMR)
library(parallel)
library(foreach)
library(doParallel)
registerDoParallel(16) # Core number

exp_pre <- fread('')#input from pQTLs datasets
exp <- exp_pre[exp_pre$cis_trans.exposure=='cis',]

outcome_path <- ''#input from full summary data for outcome phenotype
outcome_name <- 'All vs. control'#example phenotype
outcome_in <- fread(outcome_path)
colnames(outcome_in)
outcome_data <- format_data(outcome_in, 
                            type = 'outcome',
                            snps = exp$SNP,
                            snp_col = 'RsID',
                            beta_col = 'Effect',
                            se_col = 'StdErr',
                            effect_allele_col = 'Allele1',
                            other_allele_col = 'Allele2',
                            eaf_col = 'Freq1',
                            pval_col = 'Pvalue',
                            samplesize_col = 'N',
                            chr_col = 'CHR',
                            pos_col = 'POS')

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = outcome_data)
dat <- subset(dat0, dat0$mr_keep == TRUE)
length(unique(dat$gene.exposure))

mr_parallel <- function(harmonized_data)
{
  exp_list <- unique(harmonized_data$id.exposure)
  
  mr_res <- foreach (i = 1:length(exp_list)) %dopar%
    {
      MR_methods <- c('mr_wald_ratio', 
                      'mr_two_sample_ml', 
                      'mr_weighted_median', 
                      'mr_ivw_mre', 
                      'mr_ivw_fe')
      TwoSampleMR::mr(harmonized_data[harmonized_data$id.exposure == exp_list[i],],method_list = MR_methods)
    }
  return(do.call(rbind, mr_res))
}
or_res <- mr_parallel(dat)
or_res_IVWre_3 <- subset(or_res, method == 'Inverse variance weighted (multiplicative random effects)')
or_res_IVWre_3 <- or_res_IVWre_3[or_res_IVWre_3$nsnp>2,]
or_res_IVWre_2 <- subset(or_res, method == 'Inverse variance weighted (fixed effects)')
or_res_IVWre_2 <- or_res_IVWre_2[or_res_IVWre_2$nsnp<3,]
or_res_WR <- subset(or_res, method == 'Wald ratio')
or_res_all <- rbind(or_res_IVWre_3,or_res_IVWre_2,or_res_WR)

write.table(or_res, file = paste0('OR_all_', outcome_name, '_pQTL_only_cis.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
write.table(or_res_all, file = paste0('OR_all_IVWre_', outcome_name, '_pQTL_only_cis.txt'), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
