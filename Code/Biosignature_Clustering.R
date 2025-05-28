library(dplyr)         # for dataframe computation
library(MALDIquant)    # for chemometrics processing
library(caret)         # for machine learning
library(mlr3)          # for machine learning
library("mlr3verse")   # for machine learning
library("mlr3learners")# for machine learning
library("mlr3tuning")  # for machine learning
library("data.table")  # for rbindlist
library("ggplot2")     # for plots


###Data preparation
###Data preparation
setwd("~/Downloads/Data for Paleobio Project 12Sep24 //") #Reading in the data is based on the code found in: "How to import multiple .csv files simultaneously in R and create a data frame" by Dahna, A., datascience+, 03/09/2019
#https://datascienceplus.com/how-to-import-multiple-csv-files-simultaneously-in-r-and-create-a-data-frame/
species_names <- list.files()
z=lapply(species_names, read.delim,skip=4)  #read in all of 305 datasets

setwd("~/Downloads/Data for Paleobio Project 12Sep24 //") #Reading in the data is based on the code found in: "How to import multiple .csv files simultaneously in R and create a data frame" by Dahna, A., datascience+, 03/09/2019
#https://datascienceplus.com/how-to-import-multiple-csv-files-simultaneously-in-r-and-create-a-data-frame/
species_names_test <- list.files()
z2=lapply(species_names_test, read.delim,skip=4)  #read in all of 305 datasets

#LIV(L) = modern living (modern biotic)
#SYN(S) = synthetic mixtures (recent abiotic)
#RMH(R) = paleo samples (old biotic)
#CW(C) = ancient coal/wood/petroleum (old biotic)
#MET(M) = meteorite (old abiotic)

species=c(rep("C",65),rep("L",123),rep("M",42),rep("R",45),rep("S",30))
species=as.factor(species)
#species2=c(rep("B",188),rep("A", 42), rep("B", 45), rep("A", 30))
#species2=as.factor(species2)
#species3=c(rep("C",65),rep("B", 123), rep("A", 42), rep("C", 45),rep("A", 30))
#species3=as.factor(species3)


NN=700 #number of m/z values
mass=seq(50,NN,1) #m/z 
MM=3240 #number of scans
N=length(species_names)



M=list()
for(i in 1:N){
  colnames(z[[i]])="mass"
  #remove commas
  z[[i]]=data.frame(do.call("rbind", strsplit(as.character(z[[i]]$mass), ",", fixed = TRUE)))
  z[[i]]=data.frame(lapply(z[[i]],as.numeric))
  colnames(z[[i]])=c("scan",as.character(seq(50,NN,1)))
  z[[i]]=z[[i]] %>% slice(1:MM)       #Selects the first MM rows
  M[[i]]=z[[i]]
}

#saveRDS(M,"C:/Users/ghystad/Documents/Templeton_grant_research/GC_MS_R_files/Templeton_R_files/M.RDS")
#M=readRDS("C:/Users/ghystad/Documents/Templeton_grant_research/GC_MS_R_files/Templeton_R_files/M.RDS")


################
#y=species_n
#y=factor(y,labels=c("A","B"))

######################################################################
#Preprocessing
#Detect the significant peaks as local max above four times the signal to noise ratio

#number of m/z values to use
MZ=151
#Create Chromatograms for each sample and each m/z value inside each sample
Sample.list=list()
for (i in 1:N){
  S=list()
  for(j in 1:MZ){
    S[[j]] = createMassSpectrum(mass=seq(1,MM,1), intensity=unlist(M[[i]][,j+1]),metaData=list(name="Chrom"))  
  }
  
  chrom = transformIntensity(S,method="sqrt") #stabilize the variance
  chrom = smoothIntensity(chrom, method="MovingAverage",halfWindowSize=5)#Smooth the data
  chrom = removeBaseline(chrom, method="SNIP") #Remove the baseline
  
  # Put the processed chromatograms back into a dataframe
  Processed.chrom.list=list()
  for (k in 1:MZ){
    Processed.chrom.list[[k]] = as.numeric(intensity(chrom[[k]]))
  }
  
  Processed.Mass_dataframe = as.data.frame(do.call(rbind, Processed.chrom.list))
  Ma=max(Processed.Mass_dataframe)
  Mi=min(Processed.Mass_dataframe)
  #Normalize across sample
  Processed.Mass_dataframe=t((Processed.Mass_dataframe-Mi)/(Ma-Mi))
  Processed.Mass_dataframe=as.data.frame(Processed.Mass_dataframe)
  S2=list()
  for(t in 1:MZ){
    S2[[t]] = createMassSpectrum(mass=seq(1,MM,1), intensity=Processed.Mass_dataframe[,t],metaData=list(name="Chrom_normalized"))  
  }
  
  peaks = detectPeaks(S2, method="MAD", halfWindowSize=20, SNR=4)
  peak.list=list()
  for (tt in 1:MZ){
    v=numeric(MM)
    scan.number=mass(peaks[[tt]])
    v[scan.number] = intensity(peaks[[tt]])
    peak.list[[tt]] = v
  }
  Processed.peaks = t(as.data.frame(do.call(rbind, peak.list)))
  row.names(Processed.peaks)=c(paste0("R", 1:MM))
  colnames(Processed.peaks)=c(paste0("M", 50:(MZ+50-1)))
  Sample.list[[i]] = Processed.peaks
  print(i)
  print(species_names[i])
}

mass.scan.list=list()
for(i in 1:N){
  sm=as.numeric(unlist(Sample.list[[i]]))     #scan and mass spectrum for each sample
  mass.scan.list[[i]]=sm
}

#Sample vs mass/scan numbers
data.mass.scan = do.call(rbind, mass.scan.list) #Put the samples into a dataframe

Bin=as.character(seq(1,3240,1))
MS=as.character(seq(50,(MZ+50-1),1))

colnames(data.mass.scan)=paste(outer(Bin, MS, paste, sep = ';'))

#Next we use hierarchical clustering with a distance of 20 to group the scan numbers for each m/z value.
#The method is based on the paper: "Sample classification from proteing mass spectrometry by peak probability contrasts"
#by Tibshirani, R., et al. (2004), Bioinformatis, 20(17):3034-44, doi:10.1093/bioinformatics/bth357 
ml=50:(MZ+50-1)
mll=length(50:(MZ+50-1))

data.mass.scan_new=as.data.frame(ifelse(data.mass.scan>0,1,0))
MZ_name=c(paste0(";",50:(MZ+50-1)))
MZ_name2=c(paste0(50:(MZ+50-1)))
MZ_name3=c(paste0(".",50:(MZ+50-1)))

L_original=list()  # Scan numbers for each cluster per m/z value selected 
LL_dataframe=list() #The max intensities for the ith m/z value for each bin of scan numbers and each sample

for (i in 1:MZ){
  data.mass.scan_new2 = data.mass.scan_new %>% select(ends_with(MZ_name[i])) #Select columns that ends with m/z=i
  scan_name=sub(basename(colnames(data.mass.scan_new2)),pattern = MZ_name[i], replacement = "" , fixed = TRUE)
  colnames(data.mass.scan_new2) = scan_name  #Name the columns with the scan numbers
  f=function(x){
    y=sum(x)
  }
  N_elements=as.numeric(apply(as.matrix(data.mass.scan_new2),2,f)) #Count the number of  elements for each scan number across the samples
  vec=as.numeric(rep(colnames(data.mass.scan_new2), times=N_elements)) #Repeat the nonzero scan numbers with its frequency
  hc=hclust(dist(vec),method="complete")  # Hiearcial clustering on the scan numbers
  clusters=cutree(hc,h=20)                # Cluster with a distance of 20
  dataframeN=as.data.frame(data.mass.scan) %>% select(ends_with(c(paste0(MZ_name[i]))))
  L.original2=list() #Scan numbers for each cluster for the ith m/z value selected 
  L_dataframe=list() #The max intensity for the jth cluster in the ith m/z value for each sample
  for (j in 1:length(unique(clusters))){
    indd=which(clusters==unique(clusters)[j])
    L.original2[[j]]=as.data.frame(unique(vec[indd]))
    colnames(L.original2[[j]])=MZ_name2[i]
    Lt=vec[indd]
    L_mean=round(mean(Lt))   #Mean number of scan number in each cluster
    dataframeN2=dataframeN[,unique(Lt)]
    va=if(length(unique(Lt))>1){
      apply(dataframeN2,1,max) #Find the max intensity over the scan numbers in one cluster for each sample
    }else {dataframeN2}
    DD=data.frame(matrix(ncol = 1, nrow =  N)) #The max intensity for the jth cluster in the ith m/z value for each sample
    colnames(DD)= paste(outer(as.character(L_mean),MZ_name2[i], paste, sep = ';'))
    DD[,1]=va
    L_dataframe[[j]]=DD
  }
  L_original[[i]]=L.original2
  LL_dataframe[[i]]=do.call(cbind, L_dataframe) #The max intensities for the ith m/z value for each bin of scan numbers and each sample
}

data.mass.scan33=do.call(cbind, LL_dataframe) # Dataframe with the samples as rows and the intertwined m/z values and scan numbers as columns
names2=colnames(data.mass.scan33)


lg=function(x){(length(which(x>0)))/(dim(data.mass.scan33)[2])}
Perc_Nnonzero=apply(data.mass.scan33,1,lg) # Ratio of nonzero feature values for each observation

data.mass.scan33 = data.mass.scan33 %>% mutate(Perc_Nnonzero)

colnames(data.mass.scan33) = make.names(colnames(data.mass.scan33))

#Removing variables with near zero variance, nearZeroVar, and 
#correlated predictors, findCorrelation, are based on the book
#"Applied Predictive Modeling" by Kuhn, M. and Johnson, K., Springer, New York, 2013 

near.zero.variance = nearZeroVar(data.mass.scan33) #Remove features with near zero variance
data.mass.scan33 = data.mass.scan33[, -near.zero.variance] #remove variables with near zero variance
dim(data.mass.scan33) #new dimension 

#Omit correlated predictors from the samples
corr.ma=cor(data.mass.scan33)
Corr.h = findCorrelation(corr.ma, 0.85)
data.mass.scan33 <- data.mass.scan33[, -Corr.h]
dim(data.mass.scan33)  #new dimension 


#Contains all samples with the selected features and corresponding y-values (A or B)
#training_transformed=data.frame(data.mass.scan33,y=y)

clustering_transformed=data.frame(data.mass.scan33)

#________________________________

#_________
##################################

saveRDS(clustering_transformed,"Biosignature_408Samples_Processedfor_Clustering.rds")
Biosign_Training_Transformed<-readRDS("~/Downloads/Biosignature_408Samples_Processedfor_Clustering.rds")


#Biosign_Training_Transformed<-readRDS("~/Downloads/training_transformed_all305data.RDS")

#Biosign_Training_Transformed<-Biosign_Training_Transformed[,-16046]

#Biosign_Training_Transformed<-clustering_transformed

# Unsupervised Methods
# Clustering Validation and Optimal Clusters

library(cluster)
library(compiler)
gower_dist <- daisy(Biosign_Training_Transformed,
                    metric = "gower",
                    type = list())

# Calculate Silhoutte width for the fit (selecting the optimum number of clusters)
sil_width <- c(NA)

for(i in 2:15){
  pam_fit <- pam(Biosign_Training_Transformed,
                 diss = F,
                 k = i)
  sil_width[i] <- pam_fit$silinfo$avg.width
}




plot(sil_width,xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(sil_width)

library(factoextra)
library(clValid)
intern <- clValid(as.matrix(results), nClust = 2:20, 
                  clMethods = c("hierarchical","kmeans","pam"), validation = "internal",metric = "correlation",maxitems=nrow(Biosign_Training_Transformed))

intern <- clValid(Biosign_Training_Transformed, nClust = 2:20, 
                  clMethods = c("hierarchical","kmeans","pam"), validation = "internal",metric = "correlation",maxitems=nrow(Biosign_Training_Transformed))

fviz_nbclust(results, FUNcluster=pam, k.max = 20) 

fviz_nbclust(results, FUNcluster=kmeans, method="gap_stat", k.max = 20)+ theme_classic()
fviz_nbclust(results, FUNcluster=kmeans, method="silhouette", k.max = 20)+ theme_classic()

fviz_nbclust(results, FUNcluster=pam, method="gap_stat", k.max = 15)+ theme_classic()
fviz_nbclust(results, FUNcluster=pam, method="silhouette", k.max = 15)+ theme_classic()

fviz_silhouette(pam_fit) 

xxx.pca1<-prcomp(Biosign_Training_Transformed, center=FALSE, scale.=FALSE, rank. = 25)
summary(xxx.pca1) 
results <- xxx.pca1$x

# Summary
library(kableExtra)
summary(intern) %>% kable() %>% kable_styling()

pam1<-eclust(results, "pam", k=9) # factoextra::
fviz_silhouette(pam1)

rownames(results)<-Biosign_Training_Transformed$ID
kmeans1<-eclust(results, "kmeans", k=9) # factoextra::


pam_fit <- pam(results,
               diss = F,
               k = 8)

clusplot(pam_fit)

rownames(training_transformed)<-fnames2


fviz_cluster(pam_fit, results)

library(ggplot2)
library(plotly)

plot_ly(x = cluster_data$recency, 
        y = cluster_data$frequency, 
        z = cluster_data$monetary_value, 
        type = "scatter3d", 
        mode = "markers", 
        color = as.factor(cluster_data$km_cluster)) %>%
  layout(title = "",
         scene = list(xaxis = list(title = "Recency"),
                      yaxis = list(title = "Frequency"),
                      zaxis = list(title = "Monetary value")))




library(ggpubr)


#hclust_avg<-hclust()
#plot(hclust(results))

# Use hcut() which compute hclust and cut the tree
hc.cut <- hcut(results, k = 5, hc_method = "complete")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = T, rect = TRUE)
# Visualize cluster
fviz_cluster(hc.cut, ellipse.type = "convex")

rect.hclust(hclust_avg , k = 5, border = 2:8)
abline(h = 3, col = 'red')

library(caret)



#Creating confusion matrix
example <- table(as.factor(pam1$clustering), as.factor(Biosign_Training_Transformed[,16046]))

pam_fit$clustering
set.seed(123)
library(mclust)
BIC <- mclustBIC(results)
#BIC <- mclustBIC(Biosign_Training_Transformed)
plot(BIC)
summary(BIC)

mod1 <- Mclust(results, x = BIC)
summary(mod1, parameters = TRUE)

plot(mod1, what = "classification")
plot(mod1)


table(mod1$classification,as.factor(Biosign_Training_Transformed[,16046]))


Biosign_Training_Transformed$ID<-species_names
table(mod1$classification,as.factor(Biosign_Training_Transformed[,8705]))

library(dbscan)

res <- hdbscan(results, minPts = 5)

plot(res)
plot(results, col = res$cluster + 1L)
plot(res, show_flat = TRUE)


# Stacked Clusters
master.cluster <- kmeans1$cluster
slave.hierarchical <- hc.cut$cluster
slave.pam <- pam1$clustering
# Preparing the stacked clustering
stacked.clustering <- rep(NA, length(master.cluster)) 
names(stacked.clustering) <- 1:length(master.cluster)

for (cluster in unique(master.cluster)) {
  indexes = which(master.cluster == cluster, arr.ind = TRUE)
  slave1.votes <- table(slave.hierarchical[indexes])
  slave1.maxcount <- names(slave1.votes)[which.max(slave1.votes)]
   
  slave1.indexes = which(slave.hierarchical == slave1.maxcount, arr.ind = TRUE)
  slave2.votes <- table(slave.pam[indexes])
  slave2.maxcount <- names(slave2.votes)[which.max(slave2.votes)]
  
  stacked.clustering[indexes] <- slave2.maxcount
}

table(stacked.clustering,as.factor(Biosign_Training_Transformed[,16046]))


#____________________________________________

#Autoencoders

# Helper packages
library(dplyr)    # for data manipulation
library(ggplot2)  # for data visualization


# Load the library
library(keras)

# Define the encoder
encoder <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(784)) %>%
  layer_dense(units = 2, activation = 'relu')

# Define the decoder
decoder <- keras_model_sequential() %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(2)) %>%
  layer_dense(units = 784, activation = 'sigmoid')

# Connect them to create the autoencoder
autoencoder <- keras_model(inputs = encoder$input, outputs = decoder(encoder$output))


# Compile the model
autoencoder %>% compile(optimizer = 'adam', loss = 'binary_crossentropy')

# Load the data
data <- Biosign_Training_Transformed[,-8705]

# Prepare the data
x_train <- data$train$x
x_test <- data$test$x
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
x_train <- x_train / 255
x_test <- x_test / 255

# Train the model
autoencoder %>% fit(x_train, x_train, epochs = 50, batch_size = 256, validation_data = list(x_test, x_test))


# Modeling packages
library(h2o)  # for fitting autoencoders

h2o.no_progress()  # turn off progress bars
h2o.init(max_mem_size = "10g")  # initialize H2O instance

features <- as.h2o(Biosign_Training_Transformed)

# Train an autoencoder
ae1 <- h2o.deeplearning(
  x = seq_along(features),
  training_frame = features,
  autoencoder = TRUE,
  hidden = 3,
  activation = 'Tanh',
  sparse = TRUE
)

# Extract the deep features
ae1_codings <- h2o.deepfeatures(ae1, features, layer = 1)
ae1_codings


fviz_nbclust(as.matrix(ae1_codings), FUNcluster=pam, k.max = 15) 
fviz_nbclust(as.matrix(ae1_codings), FUNcluster=kmeans, method="gap_stat", k.max = 15)+ theme_classic()

eclust(as.matrix(ae1_codings), "pam", k=2)


#________________________________

library(umap)
Biosign.umap <- umap(d = Biosign_Training_Transformed)

plot.iris(Biosign.umap, as.factor(Biosign_Training_Transformed))

Biosign.umap$layout

library(splom)
library(plotly) 

Biosign.labels = Biosign_Training_Transformed
Biosign.umap = umap(Biosign_Training_Transformed, n_components = 2, random_state = 15) 
layout <- Biosign.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, Biosign_Training_Transformed) 

fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~Biosign_Training_Transformed, colors = c('green','blue','red','black','orange'), type = 'scatter', mode = 'markers')%>%  
  layout(
    plot_bgcolor = "#e5ecf6",
    legend=list(title=list(text='type')), 
    xaxis = list( 
      title = "0"),  
    yaxis = list( 
      title = "1")) 

Biosign.umap = umap(Biosign_Training_Transformed, n_components = 3, random_state = 15) 
layout <- Biosign.umap[["layout"]] 
layout <- data.frame(layout) 
final <- cbind(layout, Biosign_Training_Transformed) 

fig2 <- plot_ly(final, x = ~X1, y = ~X2, z = ~X3, color = as.factor(~Biosign_Supervised$Photosyn), colors = c('green','blue','red','black','orange')) 
fig2 <- fig2 %>% add_markers() 
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = '0'), 
                                     yaxis = list(title = '1'), 
                                     zaxis = list(title = '2'))) 

fig
fig2 

fviz_nbclust(layout, FUNcluster = fanny, k.max = 15) 
fviz_nbclust(layout, FUNcluster=kmeans, method = "gap_stat", k.max = 15)+ theme_classic()


#_______________________________

library(tsne)
library(plotly)


features <- subset(Biosign_Training_Transformed, select = -c(16046)) 

set.seed(123)
tsne <- tsne(features, initial_dims = 2)
tsne <- data.frame(tsne)
pdb <- cbind(tsne,Biosign_Training_Transformed[,16046])
options(warn = -1)
fig <-  plot_ly(data = pdb ,x =  ~X1, y = ~X2, type = 'scatter', mode = 'markers', split = ~Biosign_Supervised$Photosyn)

fig <- fig %>%
  layout(
    plot_bgcolor = "#e5ecf6"
  )

fig

#________________________________

#Metadata and features

library(readxl)
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2024-BisognatureSamples-PhotosynthesisStudy-17OCT-Prabhu-RunSamplesOnly.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-BisognatureSamples-PhotosynthesisStudy-10JAN-Prabhu-RunSamplesOnly.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-19MAR-Prabhu-RunSamplesOnly-Updated.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-25MAR-Prabhu-8categories.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-25MAR-Prabhu-9categories.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-26MAR-Prabhu-9categories-Revised.xlsx")

#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-01APR-9categories-Updated.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-01APR-PlantsVAnimals.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-02APR-LIVunknown-Updated.xlsx")
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-PhotosynthesisStudy-02APR-OnlyRMH.xlsx")
X2024_BisognatureSamples_PhotosynthesisStudy_17AUG <- read_excel("~/Downloads/2025-WongPrabhu-Model02-Samples-21MAY.xlsx")

#library(stringr)
#X2024_BisognatureSamples_PhotosynthesisStudy_17AUG$`SAMPLE #`<-str_replace(X2024_BisognatureSamples_PhotosynthesisStudy_17AUG$`SAMPLE #`,"-","")

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

# #Photosynthetic
# Biosign_Supervised_Sub <- na.omit(Biosign_Supervised[,c(3,5,21:8723)])

# #Biotic v Abiotic
# Biosign_Supervised_Sub <- na.omit(Biosign_Supervised[,c(3,7,21:8723)])
 
#Model 2 Biotic v Abiotic
Biosign_Supervised_Sub <- na.omit(Biosign_Supervised[,c(3,5:8708)])

rownames(Biosign_Supervised_Sub)
colnames(Biosign_Supervised_Sub)[1:10]

library(dplyr)
#______________________
# # First grouping
# 
# Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
#   mutate(PS_Group = case_when(
#     Photosyn %in% c(2, 3, 4, 5,10) ~ "Group A: Photosynthetic",
#     Photosyn %in% c(1,6,8,9,11) ~ "Group B: Non-photosynthetic (No synthetics)",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))
# 
# Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
#   mutate(PS_Group = case_when(
#     Photosyn %in% c(1, 2, 3) ~ "Group A: Modern Biotic",
#     Photosyn %in% c(8, 9) ~ "Group B: Abiotic",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))

Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
  mutate(PS_Group = case_when(
    Living %in% c(1) ~ "Group A: Ancient Biotic",
    Living %in% c(0) ~ "Group B: Abiotic",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))

# Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
#   mutate(PS_Group = case_when(
#     Photosyn %in% c(8,9) ~ "Group A: Abiotic",
#     Photosyn %in% c(4,5) ~ "Group B: Biotic No RMH",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))

# Biosign_Supervised_Sub <- Biosign_Supervised_Sub %>%
#   mutate(PS_Group = case_when(
#     Photosyn %in% c(1,2,3) ~ "Group A:Ancient Biotic",
#     Photosyn %in% c(8,9) ~ "Group B:Abiotic",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))


# View the result
print(Biosign_Supervised_Sub$PS_Group)

colnames(Biosign_Supervised_Sub)[1:10]

# #Run for photosyn
# Biosign_Supervised_Sub$Photosyn<-as.factor(Biosign_Supervised_Sub$PS_Group)
# summary(Biosign_Supervised_Sub$Photosyn)
# colnames(Biosign_Supervised_Sub)[1:10]
# Biosign_Supervised_Sub <- na.omit(Biosign_Supervised_Sub[,c(1:8706)])
# colnames(Biosign_Supervised_Sub)[1:10]
# Biosign_Supervised_Sub$PS_Group<-NULL

# colnames(Biosign_Supervised_Sub)[1:10]
# xxx.pca1<-prcomp(Biosign_Supervised_Sub[,-c(1,2)], center=FALSE, scale.=FALSE, rank. = 25)
# summary(xxx.pca1) 
# results <- xxx.pca1$x

# # Run for Biotic v Abiotic
# Biosign_Supervised_Sub$Living<-as.factor(Biosign_Supervised_Sub$PS_Group)
# summary(Biosign_Supervised_Sub$Living)
# colnames(Biosign_Supervised_Sub)[1:10]
# Biosign_Supervised_Sub <- na.omit(Biosign_Supervised_Sub[,c(1:8706)])
# colnames(Biosign_Supervised_Sub)[1:10]
# Biosign_Supervised_Sub$PS_Group<-NULL

# Run for Ancient Biotic v Abiotic
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

#colnames(train)

library(randomForest)
# #Run for Photosyn
# rf2 <- randomForest(train$Photosyn ~ ., mtry = 4, ntree = 1000, proximity = TRUE, data = train[,-1],localImp = T)

# # Run for Biotic
# rf2 <- randomForest(train$Living ~ ., mtry = 4, ntree = 1000, proximity = TRUE, data = train[,-1],localImp = T)

# Run for Ancient Biotic
rf2 <- randomForest(train$Living ~ ., mtry = 4, ntree = 1000, proximity = TRUE, data = train[,-1],localImp = T)

rf2

rf2$votes
write.csv(rf2$votes,"~/Downloads/Biosign_Model2_Training_Probability.csv")
# mean(rf2$votes[,1])
# mean(rf2$votes[,2])
hist(rf2$votes[,1])
# hist(rf2$votes[,2])

# For Photosyn
wrong_pred<-data.frame(File=rownames(train)[which(rf2$predicted!=train$Photosyn)], Description = train[which(rf2$predicted!=train$Photosyn),"Description"],Category="OOB",Assigned_Photosynthesis_Val=train[which(rf2$predicted!=train$Photosyn),"Photosyn"],Pred_Photosynthesis_Val=rf2$predicted[which(rf2$predicted!=train$Photosyn)])

# # For Biotic
# wrong_pred<-data.frame(File=rownames(train)[which(rf2$predicted!=train$Living)], Description = train[which(rf2$predicted!=train$Living),"Description"],Category="OOB",Assigned_Biotic_Val=train[which(rf2$predicted!=train$Living),"Living"],Pred_Biotic_Val=rf2$predicted[which(rf2$predicted!=train$Living)])

# For Photosyn
rf2_pred<-predict(rf2, test[,-1])
which(rf2_pred!=test$Photosyn)

# # For Biotic
# rf2_pred<-predict(rf2, test[,-1])
# which(rf2_pred!=test$Living)


# Test class probabilities
rf2_pred_Prob<-predict(rf2, test[,-1],type="prob")
rf2_pred_Prob
write.csv(rf2_pred_Prob,"~/Downloads/Model1_Testset_ClassProb.csv")
nrow(rf2_pred_Prob)

library(caret)
test_cm <- confusionMatrix(data=rf2_pred, reference = test$Living)
test_cm

# For Photosyn
wrong_test_pred<-data.frame(File=rownames(test)[which(rf2_pred!=test$Photosyn)], Description = test[which(rf2_pred!=test$Photosyn),"Description"],Category="Test pred",Assigned_Photosynthesis_Val=test[which(rf2_pred!=test$Photosyn),"Photosyn"],Pred_Photosynthesis_Val=rf2_pred[which(rf2_pred!=test$Photosyn)])
#rownames(test)[which(rf2_pred!=test$Photosyn)]

# # For Biotic
# wrong_test_pred<-data.frame(File=rownames(test)[which(rf2_pred!=test$Living)], Description = test[which(rf2_pred!=test$Living),"Description"],Category="Test pred",Assigned_Biotic_Val=test[which(rf2_pred!=test$Living),"Living"],Pred_Biotic_Val=rf2_pred[which(rf2_pred!=test$Living)])

#wrong_pred <-c(wrong_pred,rownames(test)[which(rf2_pred!=test$Photosyn)])
wrong_final_pred<-rbind(wrong_pred,wrong_test_pred)

write.csv(wrong_final_pred,"~/Downloads/Biosign_Model2_WronglyPred_May21.csv")
#write.csv(wrong_pred,"~/Downloads/Biosign_ModernBiotic_V_Abiotic_WronglyPred_April2.csv")
#write.csv(wrong_test_pred,"~/Downloads/Biosign_Photosyn_4v1_WronglyPred_April1.csv")

library(randomForestExplainer)
randomForestExplainer::explain_forest(rf2)

# _________________________________________________
#ROCR
library(ROCR)
library(cutpointr)

# cp <- cutpointr(test, pred, as.numeric(test$Photosyn), 
#                 method = maximize_metric, metric = sum_sens_spec)

# For training set
pred_train <- prediction(as.numeric(rf2$predicted), as.numeric(train$Living))
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

# Optimal Cutoff Point
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(roc.perf, pred))

cost.perf = performance(pred, "cost")
pred@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]

acc.perf = performance(pred, measure = "acc")
plot(acc.perf)

#______________________________________________________________

# Truly Unknown
colnames(Biosign_Supervised)

# # For photosyn
# Biosign_Supervised_Unknown <- Biosign_Supervised[,c(1,3,4,21:8723)]

#For Model2
Biosign_Supervised_Unknown <- Biosign_Supervised[,c(1,3,5,6:8708)]
#Biosign_Supervised_Unknown <- Biosign_Supervised[which(is.na(Biosign_Supervised[,c(4,21:8723)])),]

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown, "prob")

#Biosign_Supervised_Unknown$Prediction <- rf2_Unknownpred
rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

#merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Divided Photosyn")],rf2_Unknownpred)
merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Living")],rf2_Unknownpred)

write.csv(merged_unknown_results,"~/Downloads/Biosign_NoRMHBiotic_v_Abiotic_All406_UnknownPred_April7.csv")

#____________________________________________________________________

# Unknowns left after training
colnames(Biosign_Supervised)
tail(colnames(Biosign_Supervised))
colnames(Biosign_Supervised)[1:20]

tail(colnames(Biosign_Supervised_Sub))

# Biosign_Supervised_Unknown <- Biosign_Supervised %>%
#   mutate(PS_Group = case_when(
#     `Divided Photosyn` %in% c(8,9) ~ "Group A: Abiotic",
#     `Divided Photosyn` %in% c(4,5) ~ "Group B: Biotic No RMH",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))

# Biosign_Supervised_Unknown <- Biosign_Supervised %>%
#   mutate(PS_Group = case_when(
#     `Divided Photosyn` %in% c(1, 2, 3) ~ "Group A: Modern Biotic",
#     `Divided Photosyn` %in% c(8, 9) ~ "Group B: Abiotic",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))

Biosign_Supervised_Unknown <- Biosign_Supervised %>%
  mutate(PS_Group = case_when(
    Living %in% c(1) ~ "Group A: Ancient Biotic",
    Living %in% c(0) ~ "Group B: Abiotic",
    NA ~ "Other" # Optional, to catch any unexpected values
  ))

# Biosign_Supervised_Unknown <- Biosign_Supervised %>%
#   mutate(PS_Group = case_when(
#     Photosyn %in% c(2, 3, 4, 5) ~ "Group A: Photosynthetic",
#     Photosyn %in% c(1,6,8) ~ "Group B: Non-photosynthetic",
#     NA ~ "Other" # Optional, to catch any unexpected values
#   ))


colnames(Biosign_Supervised_Unknown)
Biosign_Supervised_Unknown$PS_Group

# # For Photosyn
# Biosign_Supervised_Unknown_Test <- Biosign_Supervised_Unknown[is.na(Biosign_Supervised_Unknown$PS_Group),c(1,3,5,21:8723)]

# # For Biotic
# Biosign_Supervised_Unknown_Test <- Biosign_Supervised_Unknown[is.na(Biosign_Supervised_Unknown$PS_Group),c(1,3,7,21:8723)]
# Biosign_Supervised_Unknown <- Biosign_Supervised[which(is.na(Biosign_Supervised[,c(4,21:8723)])),]

# For Model2 Biotic
Biosign_Supervised_Unknown_Test <- Biosign_Supervised_Unknown[is.na(Biosign_Supervised_Unknown$PS_Group),c(1,3,5,6:8708)]
Biosign_Supervised_Unknown <- Biosign_Supervised[which(is.na(Biosign_Supervised[,c(5,6:8708)])),]

colnames(Biosign_Supervised_Unknown_Test)

rf2_Unknownpred<-predict(rf2, Biosign_Supervised_Unknown_Test[,-c(1:3)], "prob")

#Biosign_Supervised_Unknown$Prediction <- rf2_Unknownpred
rf2_Unknownpred<-as.data.frame(rf2_Unknownpred)
rf2_Unknownpred$`Analysis File #`<-rownames(rf2_Unknownpred)

# # For Photosyn
# merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Divided Photosyn")],rf2_Unknownpred)

# For Biotic
merged_unknown_results<-merge(Biosign_Supervised_Unknown[,c("Analysis File #","Description","Living")],rf2_Unknownpred)

# # For Photosyn
# write.csv(merged_unknown_results,"~/Downloads/Biosign_Photosyn_v_NonPhotosyn_NoSyntheticsModel_185samples_UnknownPred_April10.csv")

# For Biotic
write.csv(merged_unknown_results,"~/Downloads/Biosign_Model2_109samples_UnknownPred_May21.csv")

#__________________________________________________________

#Hyperparameter tuning

# Manual Search
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
tunegrid <- expand.grid(.mtry=4)
modellist <- list()
for (ntree in c(1000, 1500, 2000, 2500)) {
  set.seed(111)
  fit <- train(train$Photosyn ~ ., data=train[,-1], method="rf", metric=metric, tuneGrid=tunegrid, trControl=control, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
  print(ntree)
}
# compare results
results <- resamples(modellist)
summary(results)
dotplot(results)

# Extend Caret
metric <- "Accuracy"
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
  predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

# train model
control <- trainControl(method="repeatedcv", number=10, repeats=3)
tunegrid <- expand.grid(.mtry=c(1:4), .ntree=c(1000, 1500, 2000, 2500))
set.seed(111)
custom <- train(train$Photosyn ~ ., data = train[,-1], method=customRF, metric=metric, tuneGrid=tunegrid, trControl=control)
summary(custom)
plot(custom)
