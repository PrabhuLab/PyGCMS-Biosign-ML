#Biosignature Model 4
#Code by Anirudh Prabhu

#Load Preprocessed data file
Biosign_Training_Transformed<-readRDS("~/Documents/Biosignature_408Samples_Processedfor_Clustering.rds")

species_names <- list.files(path = "Documents/DataforPaleobioProject12Sep24/")

library(readxl)
X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-01APR-9categories-Updated.xlsx")

Biosign_Training_Transformed$`Analysis File #` <- species_names
Biosign_Supervised<- merge(X2024_BisognatureSamples_PhotosynthesisStudy_17AUG,Biosign_Training_Transformed)

summary(Biosign_Supervised$Photosyn)
summary(Biosign_Supervised$Living)

summary(Biosign_Supervised_Sub$PS_Group)

Biosign_Supervised$Photosyn<-as.factor(Biosign_Supervised$Photosyn)
Biosign_Supervised$Photosyn<-as.factor(Biosign_Supervised$`Divided Photosyn`)
Biosign_Supervised$Living<-as.factor(Biosign_Supervised$Living)
Biosign_Supervised$Abiotic<-as.factor(Biosign_Supervised$Abiotic)
Biosign_Supervised$Fungi<-as.factor(Biosign_Supervised$Fungi)
Biosign_Supervised$Plant<-as.factor(Biosign_Supervised$Plant)
Biosign_Supervised$Microbial<-as.factor(Biosign_Supervised$Microbial)

rownames(Biosign_Supervised)<-Biosign_Supervised$`Analysis File #`
rownames(Biosign_Supervised)
colnames(Biosign_Supervised)[1:20]

#Photosynthetic
Biosign_Supervised_Sub <- na.omit(Biosign_Supervised[,c(3,5,21:8723)])

rownames(Biosign_Supervised_Sub)
colnames(Biosign_Supervised_Sub)[1:10]

library(dplyr)
#______________________

Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
  mutate(PS_Group = case_when(
    Photosyn %in% c(2, 3, 4, 5,10) ~ "Group A: Photosynthetic",
    Photosyn %in% c(1,6,8,9,11) ~ "Group B: Non-photosynthetic (No synthetics)",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))


# View the result
print(Biosign_Supervised_Sub$PS_Group)

colnames(Biosign_Supervised_Sub)[1:10]

#Run for photosyn
Biosign_Supervised_Sub$Photosyn<-as.factor(Biosign_Supervised_Sub$PS_Group)
summary(Biosign_Supervised_Sub$Photosyn)
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

#colnames(train)

library(randomForest)

rf2 <- randomForest(train$Photosyn ~ ., mtry = 4, ntree = 1000, proximity = TRUE, data = train[,-1],localImp = T)
rf2

rf2$votes
write.csv(rf2$votes,"~/Downloads/Biosign_Model4_Training_Probability.csv")
hist(rf2$votes[,1])


wrong_pred<-data.frame(File=rownames(train)[which(rf2$predicted!=train$Photosyn)], Description = train[which(rf2$predicted!=train$Photosyn),"Description"],Category="OOB",Assigned_Photosynthesis_Val=train[which(rf2$predicted!=train$Photosyn),"Photosyn"],Pred_Photosynthesis_Val=rf2$predicted[which(rf2$predicted!=train$Photosyn)])

rf2_pred<-predict(rf2, test[,-1])
which(rf2_pred!=test$Photosyn)

# Test class probabilities
rf2_pred_Prob<-predict(rf2, test[,-1],type="prob")
rf2_pred_Prob
write.csv(rf2_pred_Prob,"~/Downloads/Model1_Testset_ClassProb.csv")
nrow(rf2_pred_Prob)

library(caret)
test_cm <- confusionMatrix(data=rf2_pred, reference = test$Photosyn)
test_cm

# For Photosyn
wrong_test_pred<-data.frame(File=rownames(test)[which(rf2_pred!=test$Photosyn)], Description = test[which(rf2_pred!=test$Photosyn),"Description"],Category="Test pred",Assigned_Photosynthesis_Val=test[which(rf2_pred!=test$Photosyn),"Photosyn"],Pred_Photosynthesis_Val=rf2_pred[which(rf2_pred!=test$Photosyn)])


# # For Biotic
# wrong_test_pred<-data.frame(File=rownames(test)[which(rf2_pred!=test$Living)], Description = test[which(rf2_pred!=test$Living),"Description"],Category="Test pred",Assigned_Biotic_Val=test[which(rf2_pred!=test$Living),"Living"],Pred_Biotic_Val=rf2_pred[which(rf2_pred!=test$Living)])

#wrong_pred <-c(wrong_pred,rownames(test)[which(rf2_pred!=test$Photosyn)])
wrong_final_pred<-rbind(wrong_pred,wrong_test_pred)

write.csv(wrong_final_pred,"~/Downloads/Biosign_Model4_WronglyPred_May21.csv")

library(randomForestExplainer)
randomForestExplainer::explain_forest(rf2)

# _________________________________________________
#ROCR
library(ROCR)
library(cutpointr)

# For training set
pred_train <- prediction(as.numeric(rf2$predicted), as.numeric(train$Photosyn))
roc.perf <- performance(pred_train, measure = "tpr", x.measure = "fpr") 
plot(roc.perf, col=rainbow(10))
auc<- performance(pred_train,"auc")
auc@y.values

# For test set
predictions <- as.numeric(predict(rf2, test, type="response"))
pred <- prediction(predictions, as.numeric(test$Photosyn))
roc.perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(roc.perf, col=rainbow(10))
auc<- performance(pred,"auc")
auc@y.values

#______________________________________________________________

# All Samples
colnames(Biosign_Supervised)

# For photosyn
Biosign_Supervised_Unknown <- Biosign_Supervised[,c(1,3,4,21:8723)]

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown, "prob")

rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Divided Photosyn")],rf2_Unknownpred)
write.csv(merged_unknown_results,"~/Downloads/Biosign_NoRMHBiotic_v_Abiotic_All406_UnknownPred_April7.csv")

#____________________________________________________________________

# Unknowns left after training
colnames(Biosign_Supervised)
tail(colnames(Biosign_Supervised))
colnames(Biosign_Supervised)[1:20]

tail(colnames(Biosign_Supervised_Sub))

Biosign_Supervised_Unknown <- Biosign_Supervised %>%
  mutate(PS_Group = case_when(
    Photosyn %in% c(2, 3, 4, 5, 10) ~ "Group A: Photosynthetic",
    Photosyn %in% c(1,6,8, 9, 11) ~ "Group B: Non-photosynthetic",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))


colnames(Biosign_Supervised_Unknown)
Biosign_Supervised_Unknown$PS_Group

# For Photosyn
Biosign_Supervised_Unknown_Test <- Biosign_Supervised_Unknown[is.na(Biosign_Supervised_Unknown$PS_Group),c(1,3,5,21:8723)]

colnames(Biosign_Supervised_Unknown_Test)

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown_Test[,-c(1:3)], "prob")

rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Divided Photosyn")],rf2_Unknownpred)

# For Photosyn
write.csv(merged_unknown_results,"~/Downloads/Biosign_Photosyn_v_NonPhotosyn_147samples_UnknownPred_April10.csv")

#________________________________________________________

#Cross Validation

library(caret)
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=3)

## Train a random forest model
forest <- train(
  
  # Formula. We are using all variables to predict
  Photosyn ~ ., 
  
  # Source of data; remove the Species variable
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
