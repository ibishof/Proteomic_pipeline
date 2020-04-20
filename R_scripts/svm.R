library(e1071)

data(iris)
attach(iris)

model <- svm(diagnosis_MS ~ ., data = soma_dt)

summary(model)
print(model)
soma_validation <- soma_dt

pred <- predict(model, soma_validation)
table(pred, y)

zzzz <- as.factor(pred)
yyyyy <- ordered(zzzz)
roc_ojb_MS <- roc(response =  answers$diagnosis_MS, predictor =yyyyy )

write.csv(pred, "svm_MS.csv")
write.csv(answers, "answers-svm_ms.csv")

