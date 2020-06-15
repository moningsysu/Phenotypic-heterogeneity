library(readr)
library(dplyr)
library(parallel)
library(plyr)
library(reshape2)

setwd("~/phenotypic_heterogeneity/Supplementary/")
load("~/phenotypic_heterogeneity/Supplementary/Supplementary.Rdata")

##Fig. S1 controlling the confounding effect by mean in CV

FigureS1<-Test_trait%>%
  ggplot(aes(x=average,y=cv))+geom_point(size=0.5,alpha=0.5)+
  geom_line(aes(x=average,y=medianY),color="red",linetype=5,size=1,)+
  annotate("segment", x = 10.503398, xend = 10.503398, y =0.8858997 , yend =1.49,color="red",linetype=2,size=0.5)+
  annotate("text", x = 12, y =1.53,label="DM=0.61",color="red")+
  annotate("segment", x = 16.943945, xend = 16.943945, y=0.3546067 , yend =0.5310998,color="red",linetype=2,size=0.5)+
  annotate("text", x = 16, y =0.3,label="DM=-0.18",color="red")+
  xlab("Mean of trait value")+
  ylab("CV of trait value")


##Fig. S2  Protein-protein interaction degree


dfObs <- df70trait2Gene_NorDM %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(obsPP = mean(sort(NorDM)[51:70]))

setwd("~/phenotypic_heterogeneity/Supplementary/")
PPI_DATA<-read.csv("./PPI.csv",stringsAsFactors = F)

dfMerge<-merge(dfObs,PPI_DATA,by.x="Gene",by.y = "ORF")

cor.test(dfMerge$obsPP,dfMerge$AMS_PPI_degree,method = "s")

ci.up <- function(x){mean(x) + sd(x)/sqrt(length(x))};
ci.down <- function(x){mean(x) - sd(x)/sqrt(length(x))};

FigureS2<-ggplot(dfMerge,aes(x=type,y=obsPP)) +
  stat_summary(fun.y=mean,geom=c("point")) +
  stat_summary(fun.ymax=ci.up,fun.ymin=ci.down,geom=c("errorbar"))+
  xlab("Protein-protein interaction degree")+
  ylab("Phenotypic potential")+
  scale_x_discrete(breaks=c("A", "B", "C","D","E","F","G","H","I","J"), labels=c("1-2", "3-5","6-9","10-14","15-19","20-24","25-30","31-45","46-80","81+"))+
  annotate("text",x=2,y=1.45,parse = TRUE,label="'ρ = 0.12,' ~ italic(P) < 10^-9")




##FigureS3 

##get yeast trait number

load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")

Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)

dfCapPot2Trait_Num<-dfGene2Trait%>%
  dplyr::filter(.,Trait %in%  unique(df70trait2Gene_NorDM$Trait))%>%
  dplyr::filter(.,Gene %in% Co_Gene)%>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(numAsCap = length(which(Type == "capacitor")),
                   numAsPot = length(which(Type == "potentiator"))) 


##get Human drug number

Homology_Cap<-filter(Yeast_Human_One2One_ortholog,Gene.stable.ID %in% toupper(Capacitor$Gene))
Homology_Pot<-filter(Yeast_Human_One2One_ortholog,Gene.stable.ID %in% toupper(Potentiator$Gene))


#FigureS3 a
setwd("~/phenotypic_heterogeneity/Supplementary/Each_homologous_gene/")
Cap<-list.files("./Capacitors/")
CapDATA<-mclapply(1:length(Cap), function(k){read_csv(paste("./Capacitors/",Cap[k],sep = ""))},mc.cores = 10)
Cap_EACH<-lapply(1:length(Cap), function(k){
  Each_CAP<-CapDATA[[k]]
})%>%rbind.fill()

Sig_Cap<-filter(Cap_EACH,Pva_G<0.05)  ##Gener_Each NHNL
df_Sig_Cap<-table(Sig_Cap$Gene)%>%as.data.frame()
names(df_Sig_Cap)<-c("Human","Num_drug")

Human_Cap_Num<-Homology_Cap%>%
  dplyr::mutate(Human=paste("exp",Homology_Cap$Human.gene.name,sep = ""))%>%
  merge(.,df_Sig_Cap,by="Human")

Cap_NumberDATA<-dfCapPot2Trait_Num%>%
  dplyr::filter(.,Gene %in% Capacitor$Gene)%>%
  dplyr::mutate(Yeast=toupper(Gene))%>%
  merge(.,Human_Cap_Num,by.x="Yeast",by.y="Gene.stable.ID")%>%
  dplyr::select(.,Yeast,numAsCap,Human.gene.name,Num_drug)

cor.test(Cap_NumberDATA$numAsCap,Cap_NumberDATA$Num_drug,method = "s")


FigureS3_a<-Cap_NumberDATA%>%
  ggplot(aes(x=numAsCap,y=Num_drug,color="#F8766D"))+
  geom_point()+
  stat_smooth(method = "lm",se = F,color="#666666")+
  annotate("text",x=6,y=6150,parse = TRUE,label= "'ρ = 0.335,' ~ italic(P)< 0.03")+
  xlab("Number of traits regulated by each yeast stabilizer")+
  ylab("Number of drugs enhanced by \n each human stabilizers ")+
  scale_color_discrete(name="Gene type",labels=("#F8766D"="Capacitor"))+
  theme(legend.position="none")


#FigureS3 b
Pot<-list.files("./Potentiators/")
PotDATA<-mclapply(1:length(Pot), function(k){read_csv(paste("./Potentiators/",Pot[k],sep = ""))},mc.cores = 10)
Pot_EACH<-lapply(1:length(Pot), function(k){
  Each_POT<-PotDATA[[k]]
})%>%rbind.fill()

Sig_Pot<-filter(Pot_EACH,Pva_L<0.05)  ##Gener_Each NHNL
df_Sig_Pot<-table(Sig_Pot$Gene)%>%as.data.frame()
names(df_Sig_Pot)<-c("Human","Num_drug")

Human_Pot_Num<-Homology_Pot%>%
  dplyr::mutate(Human=paste("exp",Homology_Pot$Human.gene.name,sep = ""))%>%
  merge(.,df_Sig_Pot,by="Human")

Pot_NumberDATA<-dfCapPot2Trait_Num%>%
  dplyr::filter(.,Gene %in% Potentiator$Gene)%>%
  dplyr::mutate(Yeast=toupper(Gene))%>%
  merge(.,Human_Pot_Num,by.x="Yeast",by.y="Gene.stable.ID")%>%
  dplyr::select(.,Yeast,numAsPot,Human.gene.name,Num_drug)

cor.test(Pot_NumberDATA$numAsPot,Pot_NumberDATA$Num_drug,method = "s")

FigureS3_b<-Pot_NumberDATA%>%
  ggplot(aes(x=numAsPot,y=Num_drug,color="#00BFC4"))+
  geom_point()+
  stat_smooth(method = "lm",se = F,color="#666666")+
  annotate("text",x=6,y=4000,parse = TRUE,label= "'ρ = 0.051,' ~ italic(P)< 0.8")+
  xlab("Number of traits regulated by each yeast diversifier")+
  ylab("Number of drugs weakened by \n each human diversifier ")+
  scale_color_manual(values="#00BFC4",name="Gene type",
                     breaks="#00BFC4",labels="Potentiator")+
  scale_x_continuous(breaks = c(seq(0,16,4)))+
  theme(legend.position="none")

#FigureS3 C

##Capacitor

setwd("~/phenotypic_heterogeneity/Supplementary/Each_homologous_gene/")
Cap<-list.files("./Capacitors/")
CapDATA<-mclapply(1:length(Cap), function(k){read_csv(paste("./Capacitors/",Cap[k],sep = ""))},mc.cores = 10)
Cap_EACH<-lapply(1:length(Cap), function(k){
  Each_CAP<-CapDATA[[k]]
})%>%rbind.fill()

S3C_Cap<-Cap_EACH%>%
  dplyr::group_by(Gene)%>%
  dplyr::summarise(Num_Sensitive=length(which(Pva_G<0.05)),
                   Num_Resistant=length(which(Pva_L<0.05)),type="Capacitor")


##Potentiators
Pot<-list.files("./Potentiators/")
PotDATA<-mclapply(1:length(Pot), function(k){read_csv(paste("./Potentiators/",Pot[k],sep = ""))},mc.cores = 10)
Pot_EACH<-lapply(1:length(Pot), function(k){
  Each_POT<-PotDATA[[k]]
})%>%rbind.fill()

S3C_Pot<-Pot_EACH%>%
  dplyr::group_by(Gene)%>%
  dplyr::summarise(Num_Sensitive=length(which(Pva_G<0.05)),
                   Num_Resistant=length(which(Pva_L<0.05)),type="Potentiator")


S3C<-rbind(S3C_Cap,S3C_Pot)
cor.test(S3C$Num_Sensitive,S3C$Num_Resistant,method = "s")

FigureS3_c<-S3C_Cap%>%
  rbind(S3C_Pot)%>%
  ggplot(aes(x=Num_Resistant,y=Num_Sensitive))+
  geom_point()+
  geom_smooth(method = "lm",se=FALSE,color="red")+
  xlab("Number of drugs weakened by human diversifiers")+
  ylab("Number of drugs enhanced by  \n  human stabilizers")+
  scale_y_continuous(limits = c(0,6300),breaks = seq(0,6000,2000))+
  annotate("text", x=3000, y=6000,
           parse = TRUE,label="'ρ = -0.61,'  ~ italic(P) < ~ 10^-7")





##Draw  DATA save in Supplementary.Rdata

