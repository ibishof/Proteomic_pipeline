library(readr)
#install.packages("Rtsne")
library(Rtsne)

setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\soma_5k\\SNR\\meta_analysis")
training_set <- read.csv("midnightblue_pca_other_diagnosis_MS_only.csv")
training_set$label <- as.factor(training_set$Diagnosis_MS)
dim(training_set)

numTrain <- 50000
set.seed(3449)
#rows <- sample(1:nrow(training_set), numTrain)
#train <- training_set[rows,]


set.seed(1) # for reproducibility
tsne <- Rtsne(training_set[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 5000)

colors = rainbow(length(unique(training_set$label)))
names(colors) = unique(training_set$label)
par(mgp=c(2.5,1,0))
plot(tsne$Y, t='n', main="Immune Module tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels=training_set$label, col=colors[training_set$label])

