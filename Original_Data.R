library(dplyr)
library(parallel)
library(plyr)
library(reshape2)
library(cluster)


##1 Calculated DM of 220 traits X 4718 mutant genes
load("~/phenotypic_heterogeneity/Original_Data/220Trait_4718Gene_CV_Mean.Rdata")
df220Trait_4718Gene_DM<-df220Trait_4718Gene_CV_Mean%>%
  dplyr::group_by(Trait)%>%
  dplyr::arrange(Mean)%>%
  dplyr::mutate(Median_CV=runmed(Cv,101),
                DM=Cv-Median_CV)



##2 Using Random sampling to identify phenotypic capacitor and potentiators

##2.1 Randomly sampled data comes from wildtype
##2.2 Get Trait A in 4718 mutant genes
##2.3 Get Trait A in wildtype
##2.4 Get mutant gene a in Trait A 
##2.5 Sampling mutant gene a in wildtype, sample size is equal to the size of mutant gene a in Trait A 

load("~/phenotypic_heterogeneity/Original_Data/WildType.Rdata")
Traits<-unique(df220Trait_4718Gene_DM$Trait)
Genes<-unique(df220Trait_4718Gene_DM$Gene)
RandomSample<-mclapply(1:220,function(i){
  
  
  Each_Trait<-df220Trait_4718Gene_DM%>%
    dplyr:: filter(.,Trait %in% Traits[i])
  
  Wild_Trait<-dfWildType%>%
    dplyr:: filter(.,Trait %in% Traits[i]) 
  
  All_Gene_Samples<-lapply(1:4718, function(j){
    
    Each_Gene<-Each_Trait%>%
      dplyr:: filter(.,Gene %in% Genes[j])  
    
    A_Gene_Samples<-mclapply(1:1000,function(k){
      
      Each_Sample<-sample(Wild_Trait$value,size=Each_Gene$Cell_Number,replace = F)
      mySample<-data.frame(Trait=Traits[i],Gene=Genes[j],
                           Mean=Each_Gene$Mean,Cv=Each_Gene$Cv,
                           MedianCV=Each_Gene$Median_CV,
                           DM=Each_Gene$DM,Time=k,
                           Sample_Mean=mean(Each_Sample),
                           Sample_CV=sd(Each_Sample)/mean(Each_Sample))
      
      
    })%>%rbind.fill()
    return( A_Gene_Samples) 
    
  })%>%rbind.fill()
  
  return(All_Gene_Samples)
  
},mc.cores=10)%>%rbind.fill()

##2.6 Calculating Random_DM
lapply(1:length(Traits), function(i){
  
  A_Sample_Trait<-dplyr::filter(RandomSample,Trait%in% Traits[i])
  
  Result_A_Sample_Trait<-lapply(1:1000, function(j){
    gc()
    Sample_4718Genes<-A_Sample_Trait%>%
      dplyr::filter(.,Time %in% j)%>%
      dplyr::arrange(.,Mean)%>%
      dplyr::mutate(Sort=1:4718)
    
    dfRandomDM<-mclapply(1:4718,function(k){
      
      if(min(which(Sample_4718Genes$Sample_Mean[k]< Sample_4718Genes$Mean))!=1){
        loca<-max(which(Sample_4718Genes$Sample_Mean[k]>Sample_4718Genes$Mean)) 
        RandomCV<-Sample_4718Genes$MedianCV[which(Sample_4718Genes$Sort==loca)]
        
      }
      else {RandomCV<-Sample_4718Genes$MedianCV[which(Sample_4718Genes$Sort==1)]
      loca<-1}
      
      myfile<-data.frame(Trait=Traits[i],
                         Gene=Sample_4718Genes$Gene[k],
                         DM=Sample_4718Genes$DM[k],
                         Time=j,
                         Random_DM=Sample_4718Genes$Sample_CV[k]-RandomCV)
      
      
      
      return(myfile)
    },mc.cores = 10)%>%rbind.fill()
    
  })%>%rbind.fill()
  
  save(Result_A_Sample_Trait,file =paste( "~/phenotypic_heterogeneity/Original_Data/DM.Results/Random_DM/",Traits[i],".Rdata",sep = ""))
  
})


##2.7 To check the significance of DM 
##2.8 To identify phenotypic capacitor and potentiators
setwd("~/phenotypic_heterogeneity/Original_Data/DM.Results/Random_DM/")
myfiles<-list.files()
dfGene2Trait<-mclapply(1:220, function(i){
  load(myfiles[[i]])
  Result_A_Sample_Trait<-FII
  
  Pvalue_Data<-Result_A_Sample_Trait%>%
    dplyr::group_by(Trait,Gene,DM)%>%
    dplyr::summarise(Rank= rank(c(DM[1],Random_DM))[1] / length(Random_DM))
  
  Result_GeneType<-Pvalue_Data%>%
    dplyr::group_by(Trait,Gene)%>%
    dplyr::summarise(Type=if_else(Rank==0.001,"potentiator",
                                  ifelse(Rank >1-0.001,"capacitor","Other")))%>%
    dplyr::filter(.,Type %in% c("capacitor","potentiator"))
  
  return(Result_GeneType)                  
  
  
  
},mc.cores=10)%>%rbind.fill()



## 3 Get Trait importance (TI) in 220 traits;
## 3.1 Get 4718 mutant genes expression and Dn/Ds
## 3.2 One_to_One ortholog between yeast and human
load("~/phenotypic_heterogeneity/Original_Data/Yeast.Rdata")



## 4  Using PAM to cluster the traits and chose 70 representative traits from the clusters.
df220Trait2Gene_NorDM <-df220Trait_4718Gene_DM%>%
  dplyr::group_by(Trait)%>%
  dplyr::mutate(NorDM=(DM-mean(DM))/sd(DM))%>%
  dplyr::select(.,1,2,8)

#4.1 cluster the traits 
df220Trait2Gene_NorDM_1<-reshape2::dcast(df220Trait2Gene_NorDM,Trait ~ Gene)
Pam_Result<-pam(df220Trait2Gene_NorDM_1,70)
Traits220<-unique(df220Trait2Gene_NorDM_1$Trait)
dfCluster<-data.frame(Trait=Traits220,Trait_ID=1:220,Cluter=Pam_Result$clustering,stringsAsFactors = F)

df70trait2Gene_NorDM<-dfCluster%>%
  dplyr::filter(.,Trait_ID %in% Pam_Result$id.med)%>%
  merge(., df220Trait2Gene_NorDM,by="Trait")%>%
  dplyr::select(.,1,4,5)


##5. At a FDR of 12%, we identified potentiators and capacitors with significantly lower and higher phenotypic potentials than expected

##5.1 capacitors

dfObs <- df70trait2Gene_NorDM %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(obsPP = mean(sort(NorDM)[51:70]))

dfEachFake <- df70trait2Gene_NorDM %>%
  dplyr::group_by(Trait) %>%
  do({
    myDf <- .;
    mclapply(1:100,function(x) {
      myDf %>%
        mutate(fakeDm = sample(NorDM),
               s=x) %>%
        return()
    },mc.cores = 20) %>% rbind.fill();
  }) %>%
  group_by(s,Gene) %>%
  dplyr::summarise(fakePP = mean(sort(fakeDm)[51:70]))

dfMeanFake <- dfEachFake %>%
  dplyr::group_by(s) %>%
  dplyr::mutate(rankPP = rank(fakePP)) %>%
  dplyr::group_by(rankPP) %>%
  dplyr::summarise(fakeMeanPP = mean(fakePP));

dfRes <- dfMeanFake$fakeMeanPP %>%
  unique %>%
  lapply(function(x){
    nFalsePos = length(which(dfMeanFake$fakeMeanPP>=x));
    nPos = length(which(dfObs$obsPP>=x));
    data.frame(
      threshold = x,
      nFalsePos = nFalsePos,
      nPos = nPos,
      rateFalsePos = nFalsePos / 4718,
      nTruePos = nPos - nFalsePos) %>%
      return()
  }) %>% rbind.fill()

FDR<-min(which(dfRes$nTruePos/dfRes$nPos>=0.875))
Capacitor<-filter(dfObs,obsPP > sort(dfRes$threshold)[FDR]) 


##5.2 potentiators 

dfObs_1 <- df70trait2Gene_NorDM %>%
  group_by(Gene) %>%
  dplyr::summarise(obsPP = mean(sort(NorDM)[1:20]))

dfEachFake_1 <- df70trait2Gene_NorDM %>%
  group_by(Trait) %>%
  do({
    myDf <- .;
    lapply(1:100,function(x) {
      myDf %>%
        mutate(fakeDm = sample(NorDM),
               s=x) %>%
        return()
    }) %>% rbind.fill();
  }) %>%
  group_by(s,Gene) %>%
  dplyr::summarise(fakePP = mean(sort(fakeDm)[1:20]))

dfMeanFake_1 <- dfEachFake_1 %>%
  dplyr::group_by(s) %>%
  dplyr::mutate(rankPP = rank(fakePP)) %>%
  dplyr::group_by(rankPP) %>%
  dplyr::summarise(fakeMeanPP = mean(fakePP));

dfRes_1 <- dfMeanFake_1$fakeMeanPP %>%
  unique %>%
  lapply(function(x){
    nFalsePos = length(which(x>dfMeanFake_1$fakeMeanPP));
    nPos = length(which(x>dfObs_1$obsPP));
    data.frame(
      threshold = x,
      nFalsePos = nFalsePos,
      nPos = nPos,
      rateFalsePos = nFalsePos / 4718,
      nTruePos = nPos - nFalsePos) %>%
      return()
  }) %>% rbind.fill();

FDR_1<-max(which(dfRes_1$nTruePos/dfRes_1$nPos>=0.875))
Potentiator<-filter(dfObs_1,obsPP < sort(dfRes_1$threshold)[FDR_1]) 


##6. Save in Original.Rdata
##From 1.df220Trait_4718Gene_DM
##From 2.dfGene2Trait
##From 3.Yeast.Rdata (Contains data of Expression,DNDS ,Trait importance,and One_to_One ortholog)
##From 4.df70trait2Gene_NorDM
##From 5.Capacitor and Potentiator

save(df220Trait_4718Gene_DM,dfGene2Trait,dfGene2DndsExp,dfTrait2TI,Yeast_Human_One2One_ortholog,df70trait2Gene_NorDM,Capacitor,Potentiator,
     file = "~/phenotypic_heterogeneity/Original_Data/Original.Rdata")


load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")

