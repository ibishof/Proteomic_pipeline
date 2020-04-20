
setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\soma_5k\\SNR\\meta_analysis")
UT<-read.csv("top_feature_chi_input_claire.csv", row.names = 1)
results <- chisq.test(UT)

write.csv(results$residuals, "cliare_top_features_chi_square_results.csv")
