library(rcellminer)
library(rcellminerData)
library(dplyr)
library(parallel)
library(plyr)
library(reshape2)
library(pbapply)

##Draw  DATA save in Figure_5.Rdata

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

##get exp zscore and drug zscore
dfDrugAct <- exprs(getAct(rcellminerData::drugData)) %>% as.data.frame();
molData <- getMolDataMatrices();
expGenes <- removeMolDataType(rownames(molData[["exp"]]));
dfExpr <- molData[["exp"]][names(expGenes),] %>% as.data.frame()
rownames(dfExpr) <- expGenes;
dfExpr$gene <- expGenes;

##gene
hsStab <- Homology_Cap$Human.gene.name %>% as.character(); ## stabilizer
hsDiver <- Homology_Pot$Human.gene.name %>% as.character(); ## diversifier

hsDiver.withExpr <- intersect(hsDiver,dfExpr$gene);
dfExpr.diver <- dfExpr %>%
  filter(gene %in% hsDiver.withExpr) %>%
  melt(id.vars="gene") %>%
  group_by(gene) %>%
  plyr::rename(c("variable" = "cell", "value"="exp"));

hsStab.withExpr <- intersect(hsStab,dfExpr$gene);
dfExpr.stab <- dfExpr %>%
  filter(gene %in% hsStab.withExpr) %>%
  melt(id.vars="gene") %>%
  group_by(gene) %>%
  plyr::rename(c("variable" = "cell", "value"="exp"));

dfObs <- mclapply(
  rownames(dfDrugAct),
  mc.cores = 90,
  #mc.preschedule = F,
  function(thisDrug){
    myDrugAct <- t(dfDrugAct[thisDrug,]) %>% as.data.frame;
    colnames(myDrugAct) <- "act";
    myDrugAct$cell <- rownames(myDrugAct);
    myHiVsLow.diver <- myDrugAct %>%
      merge(dfExpr.diver,by=c("cell")) %>%
      group_by(gene) %>%
      do({
        myDf <- .;
        medExp <- myDf$exp %>% median(na.rm=T);
        data.frame(lowExp_LC50 = mean(myDf$act[myDf$exp < medExp],na.rm=T),
                   hiExp_LC50 = mean(myDf$act[myDf$exp > medExp],na.rm=T));
      });
    rawP.diver <- wilcox.test(myHiVsLow.diver$lowExp_LC50,
                              myHiVsLow.diver$hiExp_LC50,
                              paired=T,
                              alternative = "greater")$p.value;
    
    myHiVsLow.stab <- myDrugAct %>%
      merge(dfExpr.stab,by=c("cell")) %>%
      group_by(gene) %>%
      do({
        myDf <- .;
        medExp <- myDf$exp %>% median(na.rm=T);
        data.frame(lowExp_LC50 = mean(myDf$act[myDf$exp < medExp],na.rm=T),
                   hiExp_LC50 = mean(myDf$act[myDf$exp > medExp],na.rm=T));
      });
    rawP.stab <- wilcox.test(myHiVsLow.stab$lowExp_LC50,
                             myHiVsLow.stab$hiExp_LC50,
                             paired=T,
                             alternative = "less")$p.value;
    
    data.frame(drug = thisDrug, diverP = rawP.diver, stabP = rawP.stab) %>%
      return();
  }) %>% rbind.fill();

dfObs %>%
  filter(p.adjust(diverP,method="bonferroni") < 0.05) %>%
  dim() ## 11
dfObs %>%
  filter(p.adjust(stabP,method="bonferroni") < 0.05) %>%
  dim() ## 61


##
## randomly sample the same number of genes as a group, as a control
##
nRand <- 1000;
dfMock <- mclapply(
  rownames(dfDrugAct),
  mc.cores = 90,
  #mc.preschedule = F,
  function(thisDrug){
    myDrugAct <- t(dfDrugAct[thisDrug,]) %>% as.data.frame;
    colnames(myDrugAct) <- "act";
    myDrugAct$cell <- rownames(myDrugAct);
    
    myRes.mock <- lapply(1:nRand,function(thisR){
      myDiver.mock <- sample(dfExpr$gene,length(hsDiver.withExpr));
      myStab.mock <- sample(dfExpr$gene,length(hsStab.withExpr));
      
      myHiVsLow.diver <- dfExpr %>%
        filter(gene %in% myDiver.mock) %>%
        melt(id.vars="gene") %>%
        group_by(gene) %>%
        plyr::rename(c("variable" = "cell", "value"="exp")) %>%
        merge(myDrugAct,by="cell") %>%
        group_by(gene) %>%
        do({
          myDf <- .;
          medExp <- myDf$exp %>% median(na.rm=T);
          data.frame(lowExp_LC50 = mean(myDf$act[myDf$exp < medExp],na.rm=T),
                     hiExp_LC50 = mean(myDf$act[myDf$exp > medExp],na.rm=T));
        });
      rawP.diver <- wilcox.test(myHiVsLow.diver$lowExp_LC50,
                                myHiVsLow.diver$hiExp_LC50,
                                paired=T,
                                alternative = "greater")$p.value;
      
      myHiVsLow.stab <- dfExpr %>%
        filter(gene %in% myStab.mock) %>%
        melt(id.vars="gene") %>%
        group_by(gene) %>%
        plyr::rename(c("variable" = "cell", "value"="exp")) %>%
        merge(myDrugAct,by="cell") %>%
        group_by(gene) %>%
        do({
          myDf <- .;
          medExp <- myDf$exp %>% median(na.rm=T);
          data.frame(lowExp_LC50 = mean(myDf$act[myDf$exp < medExp],na.rm=T),
                     hiExp_LC50 = mean(myDf$act[myDf$exp > medExp],na.rm=T));
        });
      rawP.stab <- wilcox.test(myHiVsLow.stab$lowExp_LC50,
                               myHiVsLow.stab$hiExp_LC50,
                               paired=T,
                               alternative = "less")$p.value;
      
      data.frame(drug = thisDrug, diverP = rawP.diver, stabP = rawP.stab, rand=thisR) %>%
        return();
    }) %>% rbind.fill();
    return(myRes.mock);
  }) %>% rbind.fill();

## raw P value
dfObs %>%
  filter(stabP < 0.05) %>%
  dim() ## 8803
dfObs %>%
  filter(diverP < 0.05) %>%
  dim() ## 2707

dfMock %>%
  group_by(rand) %>%
  filter(stabP < 0.05) %>%
  dplyr::summarise(nDrugEnhanced = length(drug)) %>%
  ungroup() %>%
  dplyr::summarise(mean(nDrugEnhanced,na.rm=T),
                   max(nDrugEnhanced,na.rm=T))
## mean = 2070, max = 2161

dfMock %>%
  group_by(rand) %>%
  filter(diverP < 0.05) %>%
  dplyr::summarise(nDrugWeakened = length(drug)) %>%
  ungroup() %>%
  dplyr::summarise(mean(nDrugWeakened,na.rm=T),
                   max(nDrugWeakened,na.rm=T))
## mean = 933, max = 1012

## plot for Fig5 C/D
toPlot.5C <- dfMock %>%
  group_by(rand) %>%
  filter(stabP < 0.05) %>%
  dplyr::summarise(nDrugEnhanced = length(drug)) %>%
  ungroup() %>%
  ggplot(aes(x = nDrugEnhanced)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 150,geom="step") +
  annotate("segment", x = 8803, xend = 8803,y=0.1,color="#F8766D", yend =0.01,size=0.5,arrow=arrow())+
  annotate("text", x = 8000, y =0.12,label="Observed = 8803")+
  xlab("Number of drugs enhanced by stabilizers")+
  scale_x_continuous(limits = c(0,9000),)+
  annotate("text",x=8000,y=0.3,parse = TRUE,label= " ~ italic(P)< 10^-3")+
  ylab("Frequency")

toPlot.5D <- dfMock %>%
  group_by(rand) %>%
  filter(diverP < 0.05) %>%
  dplyr::summarise(nDrugWeakened = length(drug)) %>%
  ungroup() %>%
  ggplot(aes(x = nDrugWeakened)) +
  stat_bin(aes(y=..count../sum(..count..)),bins = 150,geom="step") +
  annotate("segment", x = 2707, xend = 2707,y=0.1,color="#00BFC4", yend =0.01,size=0.5,arrow=arrow())+
  annotate("text", x = 2600, y =0.12,label="Observed = 2707")+
  xlab("Number of drugs weakened by diversifiers")+
  scale_x_continuous(limits = c(0,3000),)+
  annotate("text",x=2600,y=0.25,parse = TRUE,label= " ~ italic(P)< 10^-3")+
  ylab("Frequency")


##Draw  DATA save in Figure_5.Rdata
save(Draw_Cap,Draw_Pot,Homology_Cap,Homology_Pot,
     dfObs,toPlot.5C,toPlot.5D,
     file = "~/phenotypic_heterogeneity/Fig_5_Results/Figure_5.Rdata")




