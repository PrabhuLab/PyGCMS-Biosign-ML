#Biosignature_Model2

# Code by Anirudh Prabhu

#Load Preprocessed data file
Biosign_Training_Transformed<-readRDS("~/Documents/Biosignature_408Samples_Processedfor_Clustering.rds")

species_names <- list.files(path = "Documents/DataforPaleobioProject12Sep24/")

library(readxl)

X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Documents/PaleobioSignature_Project/2025-WongPrabhu-Model02-Samples-21MAY.xlsx")

Biosign_Training_Transformed$`Analysis File #` <- species_names
Biosign_Supervised<- merge(X2024_BisognatureSamples_PhotosynthesisStudy_17AUG,Biosign_Training_Transformed)

summary(Biosign_Supervised$Photosyn)
summary(Biosign_Supervised$Living)

summary(Biosign_Supervised_Sub$PS_Group)

Biosign_Supervised$Photosyn<-as.factor(Biosign_Supervised$`Divided Photosyn`)
Biosign_Supervised$Living<-as.factor(Biosign_Supervised$Living)

rownames(Biosign_Supervised)<-Biosign_Supervised$`Analysis File #`
rownames(Biosign_Supervised)
colnames(Biosign_Supervised)[1:20]

Biosign_Supervised_Sub <- na.omit(Biosign_Supervised[,c(3,5:8708)])

rownames(Biosign_Supervised_Sub)
colnames(Biosign_Supervised_Sub)[1:10]

library(dplyr)
#______________________

Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
  mutate(PS_Group = case_when(
    Living %in% c(1) ~ "Group A: Ancient Biotic",
    Living %in% c(0) ~ "Group B: Abiotic",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))

# View the result
print(Biosign_Supervised_Sub$PS_Group)

colnames(Biosign_Supervised_Sub)[1:10]

Biosign_Supervised_Sub$Living<-as.factor(Biosign_Supervised_Sub$PS_Group)
summary(Biosign_Supervised_Sub$Living)
colnames(Biosign_Supervised_Sub)[1:10]
Biosign_Supervised_Sub <- na.omit(Biosign_Supervised_Sub[,c(1:8706)])
colnames(Biosign_Supervised_Sub)[1:10]
Biosign_Supervised_Sub$PS_Group<-NULL

#__________________________________

set.seed(111) # Set Seed so that same sample can be reproduced in future also
# Now Selecting 75% of data as sample from total 'n' rows of the data  
sample <- sample.int(n = nrow(Biosign_Supervised_Sub), size = floor(.75*nrow(Biosign_Supervised_Sub)), replace = F)
train <- Biosign_Supervised_Sub[sample, ]
test  <- Biosign_Supervised_Sub[-sample, ]

colnames(train)[1:10]

library(randomForest)
rf2 <- randomForest(train$Living ~ ., mtry = 4, ntree = 1000, proximity = TRUE, data = train[,-1],localImp = T)
rf2

rf2$votes
write.csv(rf2$votes,"~/Downloads/Biosign_Model2_Training_Probability.csv")
hist(rf2$votes[,1])


wrong_pred<-data.frame(File=rownames(train)[which(rf2$predicted!=train$Living)], Description = train[which(rf2$predicted!=train$Living),"Description"],Category="OOB",Assigned_Biotic_Val=train[which(rf2$predicted!=train$Living),"Living"],Pred_Biotic_Val=rf2$predicted[which(rf2$predicted!=train$Living)])
rf2_pred<-predict(rf2, test[,-1])
which(rf2_pred!=test$Living)

# Test class probabilities
rf2_pred_Prob<-predict(rf2, test[,-1],type="prob")
rf2_pred_Prob
write.csv(rf2_pred_Prob,"~/Downloads/Model1_Testset_ClassProb.csv")
nrow(rf2_pred_Prob)

library(caret)
test_cm <- confusionMatrix(data=rf2_pred, reference = test$Living)
test_cm

wrong_test_pred<-data.frame(File=rownames(test)[which(rf2_pred!=test$Living)], Description = test[which(rf2_pred!=test$Living),"Description"],Category="Test pred",Assigned_Biotic_Val=test[which(rf2_pred!=test$Living),"Living"],Pred_Biotic_Val=rf2_pred[which(rf2_pred!=test$Living)])

wrong_final_pred<-rbind(wrong_pred,wrong_test_pred)

write.csv(wrong_final_pred,"~/Downloads/Biosign_Model2_WronglyPred_May21.csv")

library(randomForestExplainer)
randomForestExplainer::explain_forest(rf2)

# _________________________________________________
#ROCR
library(ROCR)
library(cutpointr)

# For training set
pred_train <- prediction(as.numeric(rf2$predicted), as.numeric(train$Living))
roc.perf <- performance(pred_train, measure = "tpr", x.measure = "fpr") 
plot(roc.perf, col=rainbow(10))
auc<- performance(pred_train,"auc")
auc@y.values


# For test set
predictions <- as.numeric(predict(rf2, test, type="response"))
pred <- prediction(predictions, as.numeric(test$Living))
roc.perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(roc.perf, col=rainbow(10))
auc<- performance(pred,"auc")
auc@y.values

#______________________________________________________________

# Probabilities for all 273 samples
colnames(Biosign_Supervised)

Biosign_Supervised_Unknown <- Biosign_Supervised[,c(1,3,5,6:8708)]

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown, "prob")

rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Living")],rf2_Unknownpred)

write.csv(merged_unknown_results,"~/Downloads/Biosign_NoRMHBiotic_v_Abiotic_All406_UnknownPred_April7.csv")

#____________________________________________________________________

# Unknowns left after training
colnames(Biosign_Supervised)
tail(colnames(Biosign_Supervised))
colnames(Biosign_Supervised)[1:20]

tail(colnames(Biosign_Supervised_Sub))

Biosign_Supervised_Unknown <- Biosign_Supervised %>%
  mutate(PS_Group = case_when(
    Living %in% c(1) ~ "Group A: Ancient Biotic",
    Living %in% c(0) ~ "Group B: Abiotic",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))


colnames(Biosign_Supervised_Unknown)
Biosign_Supervised_Unknown$PS_Group

# For Model2 Biotic
Biosign_Supervised_Unknown_Test <- Biosign_Supervised_Unknown[is.na(Biosign_Supervised_Unknown$PS_Group),c(1,3,5,6:8708)]
Biosign_Supervised_Unknown <- Biosign_Supervised[which(is.na(Biosign_Supervised[,c(5,6:8708)])),]

colnames(Biosign_Supervised_Unknown_Test)

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown_Test[,-c(1:3)], "prob")

rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Living")],rf2_Unknownpred)
write.csv(merged_unknown_results,"~/Downloads/Biosign_Model2_109samples_UnknownPred_May21.csv")

#Cross Validation

library(caret)
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Train a random forest model
forest <- train(
  
  # Formula. We are using all variables to predict 
  Living ~ ., 
  
  # Source of data; remove the assigned variable
  data=Biosign_Supervised_Sub[,-1], 
  
  # `rf` method for random forest
  method='rf', 
  
  # Add repeated cross validation as trControl
  trControl=repeat_cv,
  
  # Accuracy to measure the performance of the model
  metric='Accuracy')

## Print out the details about the model
forest$finalModel

forest
