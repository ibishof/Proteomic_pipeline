
#This scriptremoves columsn with 50% or more zeros. Replaces zero to NA. Log2 tranforms. Then converts the na back to zero.
setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\WGCNA\\data_sets")

proteinGroups <-read.csv("test.csv", row.names = 1)

#get half min for later
half_min <- (min(proteinGroups[proteinGroups > 0]))/2

# Change zero to NA
proteinGroups[proteinGroups == 0] <- NA

# Remove rows with 50% or more NA
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
cleandat <- delete.na(proteinGroups, ((ncol(proteinGroups))/2))

#Log 2 tranform and change NA back to zero
Log2 <-log2(cleandat)
Log2[is.na(Log2)] <- 0

#Change zero to half min
Log2[Log2 == 0] <- half_min

table(is.na(Log2))






half_min <- (min(proteinGroups[proteinGroups > 0]))/2