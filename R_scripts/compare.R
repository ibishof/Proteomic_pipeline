
library(compareDF)
library(htmlTable)
install.packages("compareDF")
install.packages("htmlTable")

soma_dt <- read.csv("soma_ratios_modules_normalized_forest_input.csv")

setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\normalization_test")
soma_new <- read.csv("old_modules_max_normalized_ratios_orginal_values_full_paste.csv")

soma_og <- as.data.frame(soma_og)
soma_new <- as.data.frame(soma_new)


csoma = compare_df(soma_og, soma_new, c("SMF_ID_or_STS2_Identifier"))

print(csoma$html_output)

csoma$html_output


