library(dplyr)
library(ggplot2)
library(plyr)
library(parallel)

##Draw  DATA save in Figure_4.Rdata


#Fig4_A
##A schematic diagram showing the different predictions for the “resource allocation” model versus the “division of labor” model.


#Fig4_B
#The correlation between number of traits capacitated by a gene and "Number of traits potentiated by the gene

load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")
# Number of Capacitor and Potentiator regulating a trait in 70 Traits
Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)

dfCapPot2Trait_Num<-dfGene2Trait%>%
  dplyr::filter(.,Trait %in%  unique(df70trait2Gene_NorDM$Trait))%>%
  dplyr::filter(.,Gene %in% Co_Gene)%>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(numAsCap = length(which(Type == "capacitor")),
                   numAsPot = length(which(Type == "potentiator"))) 

cor.test(dfCapPot2Trait_Num$numAsCap,dfCapPot2Trait_Num$numAsPot,method = "s")

#draw Fig4_B
Fig4_B<-dfCapPot2Trait_Num%>%
  ggplot(aes(x=numAsPot,y=numAsCap))+
  geom_point(size=2,alpha=0.2)+
  geom_smooth(method = "lm",se=FALSE,color="red")+
  scale_x_continuous("Number of traits diversified by a gene") +
  scale_y_continuous("Number of traits buffered by a gene",limits = c(0,25))+
  annotate("text",x=15,y=20,parse = TRUE,label= "'ρ = -0.42,' ~ italic(P)<  ~10^-20")


##Fig4_D GO analysis in dual roles genes


##Fig4_C Random sampling to test the negative correlation Fig4_B 

##Fig4_C: 1. sample with regulatory relationship shuffled (70 Traits*484 Genes)
Real_Gene2Trait_Num<-dfGene2Trait%>%
  dplyr::filter(.,Trait %in%  unique(df70trait2Gene_NorDM$Trait))%>%
  dplyr::filter(.,Gene %in% Co_Gene)%>%
  merge(.,df70trait2Gene_NorDM,by=c("Trait","Gene"))%>%
  dplyr::group_by(Trait,Gene)%>%
  mutate(New_TYPE=ifelse(is.na(Type)==TRUE,NA,
                         ifelse(nchar(Type)>9,-1,1)))%>%
  dplyr::select(Gene,Trait,New_TYPE)%>%
  reshape2::dcast(Trait ~ Gene)

Real_Gene2Trait_Num[is.na(Real_Gene2Trait_Num)]<-0
#New_Type=-1,Gene Type is potentiator
#New_Type=1, Gene Type is capacitor
#New_Type=0, Gene Type is others

SampleDATA<-Real_Gene2Trait_Num%>%
  dplyr::select(.,-Trait)%>%as.matrix()

nRow <- nrow(SampleDATA)
nCol <- ncol(SampleDATA)

Sam_Correlation=c()

for( j in c(1:1000)){
  
  Sam<-50000+10000*j;
  
  for(i in c(1:Sam)) {
    pickRow <- sample.int(nRow,size=2);
    pickCol <- sample.int(nCol,size=2);
    SampleDATA[pickRow,pickCol] <- SampleDATA[rev(pickRow),rev(pickCol)]; ## randomly upset once
  } 
  
  Test_Sam<-SampleDATA%>%as.data.frame()
  Test_Sam$Trait<-Real_Gene2Trait_Num$Trait
  MyDATA<-Test_Sam%>%
    reshape2::melt(id.var="Trait")%>%
    dplyr::filter(value %in% c(1,-1))%>%
    dplyr::filter(variable %in% Co_Gene)%>%
    dplyr::group_by(variable)%>%
    dplyr::summarise(numAsCap=length(which(value==1)),
                     numAsPot=length(which(value==-1)))
  
  Correlation<-cor.test(MyDATA$numAsCap,MyDATA$numAsPot,method = "s")
  myfile<-data.frame(Time=j,Rho=Correlation$estimate,Pvalue=Correlation$p.value)
  Sam_Correlation=rbind(Sam_Correlation,myfile)
  
}


## 2. sample with Sample with raw DM shuffled (220 Traits*4718Genes)
## reference Original_Data .R

##2.1 Get 1000 random shuffled DM 

Traits<-unique(df220Trait_4718Gene_DM$Trait)
Genes<-unique(df220Trait_4718Gene_DM$Gene)

Rough_DM<-lapply(1:1000, function(i){
  gc()
  One_Sample<-mclapply(1:220, function(k){
    
    Each_Trait<-df220Trait_4718Gene_DM%>%
      dplyr::filter(.,Trait %in% Traits[k])%>%
      dplyr::arrange(.,Mean)%>%
      dplyr::mutate(ID=1:length(Gene))
    
    
    Each_Trait_Sam<-Each_Trait%>%
      dplyr::select(.,Mean,Cv,Median_CV,DM)%>%
      dplyr::mutate(Rough_ID=sample(Each_Trait$ID,size = 4718,replace = F),
                    Rough_Time=i)
    
    #Randomly disrupt the connect between Gene and Mean/CV/MedianCV/DM.
    mydata<-Each_Trait%>%
      dplyr::select(.,Gene,ID,Cell_Number)%>%
      merge(.,Each_Trait_Sam,by.x="ID",by.y="Rough_ID")%>%
      dplyr::rename("Trait"="Trait.x",
                    "RoughMean"="Mean",
                    "RoughMedianCV"="Median_CV",
                    "RoughDM"="DM")%>%
      dplyr::select(.,Trait,Gene,RoughMean,RoughMedianCV,RoughDM,Cell_Number,Rough_Time)
    
    return(mydata)
  },mc.cores = 30)%>%rbind.fill()
  
})

##2.2 For each random DM of 220 traits X 4718 genes, the data of random sampling comes from wild type
load("~/phenotypic_heterogeneity/Original_Data/WildType.Rdata")
setwd("~/phenotypic_heterogeneity/Fig_4.Results/")

##2.3 sample from wild type 
##2.4 Calculating Random_DM
##2.5 To check the significance of DM 
##2.6 To identify phenotypic capacitor and potentiators
lapply(1:1000, function(q){
  
  Each_Rough_DM<-Rough_DM[[q]]
  gc()
  lapply(1:220,function(i){
    
    Each_Trait<-Each_Rough_DM%>%dplyr:: filter(.,Trait %in% Traits[i])
    Wild_Trait<-dfWildType%>%dplyr:: filter(.,Trait %in% Traits[i]) 
    
    Rough_DM_Sam_WT_A_Trait<-lapply(1:1000, function(j){
      
      ##2.3 sample from wild type 
      A_Gene_Samples<-mclapply(1:4718,function(k){
        
        Each_Gene<-Each_Trait%>%dplyr:: filter(.,Gene %in% Genes[k])  
        Each_Sample_Cell_Value<-sample(Wild_Trait$value,size=Each_Gene$Cell_Number,replace = F)
        
        mySample<-data.frame(Gene=Genes[k],
                             Rough_Mean=Each_Gene$RoughMean,
                             Rough_MedianCV=Each_Gene$RoughMedianCV,
                             Rough_DM=Each_Gene$RoughDM,
                             Sam_WT_Mean=mean(Each_Sample_Cell_Value),
                             Sam_WT_CV=sd(Each_Sample_Cell_Value)/mean(Each_Sample_Cell_Value))
        
        return(mySample)
        
      },mc.cores = 20)%>%rbind.fill()
      
      Sample_False_DM<-A_Gene_Samples%>%
        dplyr::arrange(.,Rough_Mean)%>%
        dplyr::mutate(Sort=1:4718)
      
      ##2.4 Calculating Random_DM
      Sample_False_DM_1<-mclapply(1:4718,function(a){
        
        if(min(which(Sample_False_DM$Sam_WT_Mean[a]< Sample_False_DM$Rough_Mean))!=1){
          loca<-max(which(Sample_False_DM$Sam_WT_Mean[a]>Sample_False_DM$Rough_Mean)) 
          RandomCV<-Sample_False_DM$Rough_MedianCV[which(Sample_False_DM$Sort==loca)]
          
        }
        else {RandomCV<-Sample_False_DM$Rough_MedianCV[which(Sample_False_DM$Sort==1)]
        loca<-1}
        
        myfile<-data.frame(Trait=Traits[i],
                           Gene=Sample_False_DM$Gene[a],
                           RoughDM=Sample_False_DM$Rough_DM[a],
                           Sam_False_Time=j,
                           False_DM=Sample_False_DM$Sam_WT_CV[a]-RandomCV)
        
        return(myfile)
      },mc.cores = 10)%>%rbind.fill()
      
      return(Sample_False_DM_1)
      
      
    })%>%rbind.fill()
    
    ## One random sampling  in a trait has been completed
    
    ###2.5 To check the significance of DM 
    Pvalue_Data<-Rough_DM_Sam_WT_A_Trait%>%
      dplyr::group_by(Trait,Gene,RoughDM)%>%
      dplyr::summarise(Rank= rank(c(RoughDM[1],False_DM))[1] / length(False_DM))
    
    ##2.6 To identify phenotypic capacitor and potentiators
    Sam_Result_GeneType<-Pvalue_Data%>%
      dplyr::group_by(Trait,Gene,RoughDM)%>%
      dplyr::summarise(Type=if_else(Rank==0.001,"potentiator",
                                    ifelse(Rank >1-0.001,"capacitor","Other")))%>%
      dplyr::filter(.,Type %in% c("capacitor","potentiator"))
    
    save(Sam_Result_GeneType,file =paste("./Sample with raw DM shuffled/","Random_",q,"-",Traits[i],".Rdata",sep = ""))
    
  })
  
})

##2.6.1 Data processing： get one random sample significant Gene type
setwd("/mnt/data/home/moning/phenotypic_heterogeneity/Fig_4.Results/Sample with raw DM shuffled/")
Sam_dfGene2Traits<-mclapply(1:1000, function(k){
  
  myfiles<-list.files(pattern =paste("Random_",k,sep = "")) # get each sample;
  
  One_SigDA<-mclapply(1:220, function(t){
    load(myfiles[t]) #get Sam_Result_GeneType
    
    One_Trait<-Sam_Result_GeneType%>%
      dplyr::select(.,Trait,Gene,Type)
    
    return(One_Trait)
  },mc.cores = 10)%>%rbind.fill()
  
},mc.cores = 10)


## 2.7 Using PAM to cluster the traits and chose 70 representative traits from the clusters.
##.2.8 At a FDR of 12%, identified potentiators and capacitors with significantly lower and higher phenotypic potentials than expected
setwd("~/phenotypic_heterogeneity/Fig_4.Results/Sample with raw DM shuffled/")

Sam_Correlation_1<-mclapply(1:1000, function(w){
  
  ##use random dm data  to cluster the traits with PAM 
  Each_Rough_DM<-Rough_DM[[q]]
  
  #get Sam_Result and Normalized RoughDM
  One_Sample_NorDM<-Each_Rough_DM%>%
    dplyr::mutate(NorDM=(RoughDM-mean(RoughDM))/sd(RoughDM))%>%
    dplyr::select(.,Trait,Gene,NorDM)
  
  df_One_Sample_NorDM<-reshape2::dcast(One_Sample_NorDM,Trait ~ Gene)
  Pam_Result<-pam( df_One_Sample_NorDM,70)
  Traits220<-unique( df_One_Sample_NorDM$Trait)
  dfCluster<-data.frame(Trait=Traits220,Trait_ID=1:220,Cluter=Pam_Result$clustering,stringsAsFactors = F)
  
  Sam70trait2Gene_NorDM<-dfCluster%>%
    dplyr::filter(.,Trait_ID %in% Pam_Result$id.med)%>%
    merge(.,One_Sample_NorDM,by="Trait")%>%
    dplyr::select(.,Trait,Gene,NorDM)
  
  ##2.8.1 capacitors
  ##2.8.2 potentiators 
  {
    dfObs <- Sam70trait2Gene_NorDM %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(obsPP = mean(sort(NorDM)[51:70]))
    
    dfEachFake <- Sam70trait2Gene_NorDM %>%
      dplyr::group_by(Trait) %>%
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
    
    
    ##2.8.2 potentiators 
    
    dfObs_1 <- Sam70trait2Gene_NorDM %>%
      group_by(Gene) %>%
      dplyr::summarise(obsPP = mean(sort(NorDM)[1:20]))
    
    dfEachFake_1 <- Sam70trait2Gene_NorDM %>%
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
    }
  
  ##Test correlation
  Sam_Co_Gene<-union(Capacitor$Gene,Potentiator$Gene) ##Sam_Co_Gene many be 0
  
  ## load Data :One random sample significant Gene type
  Each_Sam_dfGene2Traits<-Sam_dfGene2Traits[[w]]
  
  One_Random_Rho<-Each_Sam_dfGene2Traits%>%
    dplyr::filter(.,Trait %in%  unique(Sam70trait2Gene_NorDM$Trait))%>%
    dplyr::filter(.,Gene %in% Sam_Co_Gene)%>%  
    dplyr::group_by(Gene) %>%
    dplyr::summarise(numAsCap = length(which(Type == "capacitor")),
                     numAsPot = length(which(Type == "potentiator"))) 
  
  Correlation_1<-cor.test(One_Random_Rho$numAsCap,One_Random_Rho$numAsPot,method = "s")
  
  LastDATA<-data.frame(Rho=Correlation_1$estimate,
                       Pvalue=Correlation_1$p.value,
                       SamTime=w)
  
  return(LastDATA)
  
},mc.cores = 10)%>%rbind.fill()       



##draw Fig4_C


Sam_Correlation$type<-"Samples with regulatory relationship shuffled"
Sam_Correlation_1$type<-"Samples with raw DM shuffled"
Fig4_C<-Sam_Correlation%>%
  rbind(.,Sam_Correlation_1)%>%
  ggplot(aes(x=Rho,color=type))+
  scale_x_continuous(limits = c(-0.45,0.2),breaks = c(seq(-0.4,0.2,0.1)))+
  stat_bin(aes(y=..count../sum(..count..)),bins = 100,geom="step")+
  annotate("segment", x = -.42, xend = -.42, y =0.015,color="red", yend =0.002,size=0.5,arrow=arrow())+
  annotate("text", x = -.35, y =0.018,label="Observed= -0.42")+
  scale_color_manual(name=" ",
                     labels=c("Samples with raw DMs shuffled",
                              "Samples with regulatory relationships shuffled"),
                     
                     values=c("black","#FF9933"))+
  xlab("Spearman's correlation coefficient")+
  ylab("Frequency")


##Draw  DATA save in Figure_4.Rdata
save(dfCapPot2Trait_Num,Sam_Correlation,Sam_Correlation_1,file = "~/phenotypic_heterogeneity/Fig_4_Results/Figure_4.Rdata")



