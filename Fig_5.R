library(rcellminer)
library(rcellminerData)
library(dplyr)
library(parallel)
library(plyr)
library(reshape2)
library(pbapply)

##Fig5_AB Examples
setwd("~/phenotypic_heterogeneity/")
load("./Original_Data/Original.Rdata")

##get exp zscore and drug zsocre
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()
expGenes <- removeMolDataType(rownames(molData[["exp"]]))

##gene
Homology_Cap<-filter(Yeast_Human_One2One_ortholog,Gene.stable.ID %in% toupper(Capacitor$Gene))
Homology_Pot<-filter(Yeast_Human_One2One_ortholog,Gene.stable.ID %in% toupper(Potentiator$Gene))

COGENE<-union(Homology_Cap$Human.gene.name,Homology_Pot$Human.gene.name)

CO_POT<-intersect(Homology_Pot$Human.gene.name%>%as.character(),expGenes)
POTLabels <- paste0("exp", CO_POT)
POTEXP_28 <- molData[["exp"]][POTLabels,]%>%as.data.frame()
POTEXP_28$Gene<-row.names(POTEXP_28)

CO_CAP<-intersect(Homology_Cap$Human.gene.name%>%as.character(),expGenes)
CAPLabels <- paste0("exp", CO_CAP)
CAPEXP_42 <- molData[["exp"]][CAPLabels,]%>%as.data.frame()
CAPEXP_42$Gene<-row.names(CAPEXP_42)


#Fig5_A  Stabilizer
k=3780
drugnsc<-row.names(drugAct)
Onedrug<-drugAct[drugnsc[k],]%>%as.data.frame()
names(Onedrug)<-drugnsc[k]
Onedrug$Cellline<-row.names(Onedrug)

Gene42DATA<-pblapply(1:42, function(i){
  OneGe<-CAPEXP_42[i,]
  OneGene<-melt(OneGe,id.var="Gene")
  OneGene<-na.omit(OneGene)
  OneGene$ExpType <- ifelse(OneGene$value <  quantile(OneGene$value,p=c(0.5)), "low", 
                            ifelse(OneGene$value > quantile(OneGene$value,p=c(0.5)), "high", NA))
  HIGH_Cell<-filter(OneGene,ExpType=="high")
  LOW_Cell<-filter(OneGene,ExpType=="low")
  H_act<-filter(Onedrug,Cellline %in% HIGH_Cell$variable)
  L_act<-filter(Onedrug,Cellline %in% LOW_Cell$variable)
  HAC<-mean(na.omit(H_act[,1]))
  LAC<-mean(na.omit(L_act[,1]))
  SeH<-sd(na.omit(H_act[,1]))/sqrt(nrow(H_act))
  SeL<-sd(na.omit(L_act[,1]))/sqrt(nrow(L_act))
  new<-data.frame(Gene=CAPEXP_42$Gene[i],drug=drugnsc[k],
                  Mean_H=HAC,SE_H=SeH,Mean_L=LAC,SE_L=SeL)
  return(new)
},cl=20)%>%rbind.fill()

#remove dual roles genes
AGene42DATA<-Gene42DATA%>%
  dplyr::filter(.,!(Gene %in% c( "expHDLBP" ,"expTP53RK")))

wilcox.test(AGene42DATA$Mean_H,AGene42DATA$Mean_L,alternative = "greater",paired = T)


Draw_A1<-AGene42DATA%>%
  select(.,1,2,3,4)%>%
  dplyr::rename(.,"Aver"="Mean_H","SE"="SE_H")%>%
  arrange(.,Gene)%>%
  mutate(type="High expression",GID=1:40)%>%
  arrange(.,Aver)%>%
  mutate(AID=40:1,NID=(AID%/%6)+1)

Draw_A2<-AGene42DATA%>%
  select(.,1,2,5,6)%>%
  dplyr::rename(.,"Aver"="Mean_L","SE"="SE_L")%>%
  arrange(.,Gene)%>%
  mutate(type="Low  expression",GID=1:40,)

Draw_A3<-Draw_A1%>%
  select(.,1,6,7,8)%>%
  merge(.,Draw_A2,by="Gene")%>%
  select(.,1,5,6,7,8,9,3,4)%>%
  dplyr::rename(.,"GID"="GID.y")

Draw_Cap<-Draw_A1%>%
  rbind(.,Draw_A3)%>%
  dplyr::group_by(Gene,type)%>%
  dplyr::mutate(Y_Max=Aver+SE,Y_Min=Aver-SE)


Fig5_A<-Draw_Cap%>%
  ggplot(aes(x=AID,y=Aver,color=type))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=0.4)+
  geom_point(size=2)+
  scale_color_manual(name="Type of gene expression",
                     breaks = c("High expression","Low  expression"),
                     values=c("#FF9933","grey"))+
  xlab("Stabilizer")+
  ylab("Cytotoxicity of TIC10")+
  annotate("text",x=12,y=0.6,parse = TRUE,label= "'Wilcoxon signed-rank test,' ~ italic(P)< ~10^-10")


#Fig5_B  Diversifier
k=3978
drugnsc<-row.names(drugAct)
Onedrug<-drugAct[drugnsc[k],]%>%as.data.frame()
names(Onedrug)<-drugnsc[k]
Onedrug$Cellline<-row.names(Onedrug)

Gene28DATA<-pblapply(1:28, function(i){
  OneGe<-POTEXP_28[i,]
  OneGene<-melt(OneGe,id.var="Gene")
  OneGene<-na.omit(OneGene)
  OneGene$ExpType <- ifelse(OneGene$value <  quantile(OneGene$value,p=c(0.5)), "low", 
                            ifelse(OneGene$value > quantile(OneGene$value,p=c(0.5)), "high", NA))
  HIGH_Cell<-filter(OneGene,ExpType=="high")
  LOW_Cell<-filter(OneGene,ExpType=="low")
  H_act<-filter(Onedrug,Cellline %in% HIGH_Cell$variable)
  L_act<-filter(Onedrug,Cellline %in% LOW_Cell$variable)
  HAC<-mean(na.omit(H_act[,1]))
  LAC<-mean(na.omit(L_act[,1]))
  SeH<-sd(na.omit(H_act[,1]))/sqrt(nrow(H_act))
  SeL<-sd(na.omit(L_act[,1]))/sqrt(nrow(L_act))
  new<-data.frame(Gene=POTEXP_28$Gene[i],drug=drugnsc[k],
                  Mean_H=HAC,SE_H=SeH,Mean_L=LAC,SE_L=SeL)
  new$D<-new$Mean_H-new$Mean_L
  return(new)
},cl=20)%>%rbind.fill()

#remove dual roles genes
AGene28DATA<-Gene28DATA%>%
  dplyr::filter(.,!(Gene %in% c("expHDLBP" ,"expTP53RK")))

wilcox.test(Gene28DATA$Mean_H,Gene28DATA$Mean_L,alternative = "less",paired = T)

Draw_B1<-AGene28DATA%>%
  select(.,1,2,3,4)%>%
  dplyr::rename(.,"Aver"="Mean_H","SE"="SE_H")%>%
  arrange(.,Gene)%>%
  mutate(type="High expression",GID=1:26)


Draw_B2<-AGene28DATA%>%
  select(.,1,2,5,6)%>%
  dplyr::rename(.,"Aver"="Mean_L","SE"="SE_L")%>%
  arrange(.,Gene)%>%
  mutate(type="Low  expression",GID=1:26)%>%
  arrange(.,Aver)%>%
  mutate(AID=26:1)

Draw_B3<-Draw_B2%>%
  select(.,1,6,7)%>%
  merge(.,Draw_B1,by="Gene")%>%
  select(.,1,4,5,6,7,8,3)%>%
  dplyr::rename(.,"GID"="GID.y")


Draw_Pot<-Draw_B2%>%
  rbind(.,Draw_B3)%>%
  dplyr::group_by(Gene,type)%>%
  dplyr::mutate(Y_Max=Aver+SE,Y_Min=Aver-SE)

Fig5_B<-Draw_Pot%>%
  ggplot(aes(x=AID,y=Aver,color=type))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=0.4)+
  geom_point(size=2)+
  scale_color_manual(name="Type of gene expression",
                     breaks = c("High expression","Low  expression"),
                     values=c("#FF9933","grey"))+
  xlab("Diversifier")+
  ylab("Cytotoxicity of Coumarin 340")+
  scale_x_continuous(breaks = c(seq(0,25,5)))+
  annotate("text",x=8,y=0.68,parse = TRUE,label= "'Wilcoxon signed-rank test,' ~ italic(P)< ~10^-6")





##Fig5_CD

##Real value
load("~/phenotypic_heterogeneity/Fig_5_Results/Fig5CD.Rdata")

TCap<-Capacitors%>%
  dplyr::mutate(Real_P=p.adjust(Pva_G,method = "BH"))%>%
  dplyr::filter(.,Real_P < 0.05)

nrow(TCap)/42


TGen<-Potentiators%>%
  dplyr::mutate(Real_P=p.adjust(Pva_L,method = "BH"))%>%
  dplyr::filter(.,Real_P < 0.05)

nrow(TGen)/28




###Sample from all human genes

library(rcellminer)
library(rcellminerData)
library(dplyr)
library(parallel)
library(plyr)
library(reshape2)

##get exp zscore and drug zsocre
drugAct <- exprs(getAct(rcellminerData::drugData))
molData <- getMolDataMatrices()
expGenes <- removeMolDataType(rownames(molData[["exp"]]))
lapply(1:1000, function(P){
  
  HOMOCap<-sample(expGenes,42,replace = F)
  HOMOPot<-sample(expGenes,28,replace = F)
  
  PotLabels <- paste0("exp", HOMOPot)
  PotEXP_28 <- molData[["exp"]][PotLabels,]%>%as.data.frame()
  PotEXP_28$Gene<-row.names(PotEXP_28)
  
  CapLabels <- paste0("exp", HOMOCap)
  CapEXP_42 <- molData[["exp"]][CapLabels,]%>%as.data.frame()
  CapEXP_42$Gene<-row.names(CapEXP_42)
  
  ## Each random sample
  Sam_Potentiator<-mclapply(1:21121, function(k){
    gc()
    drugnsc<-row.names(drugAct)
    Onedrug<-drugAct[drugnsc[k],]%>%as.data.frame()
    names(Onedrug)<-drugnsc[k]
    Onedrug$Cellline<-row.names(Onedrug)
    
    
    Gene28DATA<-lapply(1:28, function(i){
      OneGe<-PotEXP_28[i,]
      OneGene<-melt(OneGe,id.var="Gene")
      OneGene<-na.omit(OneGene)
      OneGene$ExpType <- ifelse(OneGene$value <  quantile(OneGene$value,p=c(0.5)), "low", 
                                ifelse(OneGene$value > quantile(OneGene$value,p=c(0.5)), "high", NA))
      HIGH_Cell<-filter(OneGene,ExpType=="high")
      LOW_Cell<-filter(OneGene,ExpType=="low")
      H_act<-filter(Onedrug,Cellline %in% HIGH_Cell$variable)
      L_act<-filter(Onedrug,Cellline %in% LOW_Cell$variable)
      HAC<-mean(na.omit(H_act[,1]))
      LAC<-mean(na.omit(L_act[,1]))
      new<-data.frame(Gene=PotEXP_28$Gene[i],drug=drugnsc[k],H_ACT=HAC,L_ACT=LAC)
      return(new)
    })%>%rbind.fill()
    
    Test_G<-wilcox.test(Gene28DATA$H_ACT,Gene28DATA$L_ACT,alternative = "greater",paired = T)
    Test_L<-wilcox.test(Gene28DATA$H_ACT,Gene28DATA$L_ACT,alternative = "less",paired = T)
    Gene28DATA$Pva_G=Test_G$p.value
    Gene28DATA$Pva_L=Test_L$p.value
    
    return(Gene28DATA)
    
  },mc.cores = 10)%>%rbind.fill()
  
  Sam_Capacitor<-mclapply(1:21121, function(k){
    gc()
    drugnsc<-row.names(drugAct)
    Onedrug<-drugAct[drugnsc[k],]%>%as.data.frame()
    names(Onedrug)<-drugnsc[k]
    Onedrug$Cellline<-row.names(Onedrug)
    
    
    Gene42DATA<-lapply(1:42, function(i){
      OneGe<-CapEXP_42[i,]
      OneGene<-melt(OneGe,id.var="Gene")
      OneGene<-na.omit(OneGene)
      OneGene$ExpType <- ifelse(OneGene$value <  quantile(OneGene$value,p=c(0.5)), "low", 
                                ifelse(OneGene$value > quantile(OneGene$value,p=c(0.5)), "high", NA))
      HIGH_Cell<-filter(OneGene,ExpType=="high")
      LOW_Cell<-filter(OneGene,ExpType=="low")
      H_act<-filter(Onedrug,Cellline %in% HIGH_Cell$variable)
      L_act<-filter(Onedrug,Cellline %in% LOW_Cell$variable)
      HAC<-mean(na.omit(H_act[,1]))
      LAC<-mean(na.omit(L_act[,1]))
      new<-data.frame(Gene=CapEXP_42$Gene[i],drug=drugnsc[k],H_ACT=HAC,L_ACT=LAC)
      return(new)
    })%>%rbind.fill()
    
    
    Test_G<-wilcox.test(Gene42DATA$H_ACT,Gene42DATA$L_ACT,alternative = "greater",paired = T)
    Test_L<-wilcox.test(Gene42DATA$H_ACT,Gene42DATA$L_ACT,alternative = "less",paired = T)
    Gene42DATA$Pva_G=Test_G$p.value
    Gene42DATA$Pva_L=Test_L$p.value
    
    return(Gene42DATA)
    
  },mc.cores = 10)%>%rbind.fill()
  
  save(Sam_Capacitor,Sam_Potentiator,file =paste("~/phenotypic_heterogeneity/Fig_5_Results/Pair_whole_Gene_in ALL/","Sample_",P,"_NCI60.Rdata",sep = ""))
  
})





##drwa Figure 5 CD
##load data :sample from all human genes
setwd("~/phenotypic_heterogeneity/Fig_5_Results/Pair_whole_Gene_in ALL/")
Allfile<-list.files("./")

Sam_ALL<-mclapply(1:length(Allfile), function(k){
  
  load(Allfile[k])
  
  Num_Cap<-Sam_Capacitor%>%
    dplyr::mutate(NC=p.adjust(Pva_G,method = "BH"))%>%
    dplyr::filter(.,NC<0.05)%>%
    dplyr::summarise(Time=k,Number=length(NC)/42,type="Cap")
  
  Num_Pot<-Sam_Potentiator%>% 
    dplyr::mutate(NP=p.adjust(Pva_L,method = "BH"))%>%
    dplyr::filter(.,NP<0.05)%>%
    dplyr::summarise(Time=k,Number=length(NP)/28,type="Pot")
  
  new<-rbind(Num_Cap,Num_Pot)
  
},mc.cores = 20)%>%rbind.fill()


Fig5_C<-Sam_ALL%>%
  filter(.,type=="Cap")%>%
  ggplot(aes(x=,Number))+
  scale_x_continuous(limits = c(0,6300),breaks = seq(0,6300,1000))+
  stat_bin(aes(y=..count../sum(..count..)),bins = 30,geom="step")+
  annotate("segment", x = 6238, xend = 6238, y =0.2,color="#F8766D", yend =0.005,size=0.5,arrow=arrow())+
  annotate("text", x =5800, y =0.24,label="Observed= 6238")+
  xlab("Number of drugs enhanced by stabilizer")+
  annotate("text",x=1000,y=0.9,parse = TRUE,label= " ~ italic(P) < 10^-3")+
  ylab("Frequency")


Fig5_D<-Sam_ALL%>%
  filter(.,type=="Pot" )%>%
  ggplot(aes(x=,Number))+
  scale_x_continuous(limits = c(0,840),breaks = seq(0,840,100))+
  stat_bin(aes(y=..count../sum(..count..)),bins = 30,geom="step")+
  annotate("segment", x = 823, xend = 823,y=0.2,color="#00BFC4", yend =0.005,size=0.5,arrow=arrow())+
  annotate("text", x =780, y =0.24,label="Observed= 823")+
  xlab("Number of drugs weakened by diversifier" )+
  annotate("text",x=120,y=0.95,parse = TRUE,label= " ~ italic(P)< 10^-3")+
  ylab("Frequency")


##Draw  DATA save in Figure_5.Rdata


