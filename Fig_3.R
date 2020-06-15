library(ggsignif)
library(dplyr)
library(ggplot2)
library(plyr)
library(parallel)


##(Fig3_A,Fig3_B, P < xxx, one tailed Wilcoxon rank-sum test)

#load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")

##dN/dS
setwd("~/phenotypic_heterogeneity/")
dfGene2dNdS<-read.csv("./Fig_3_Results/dNdS.csv",stringsAsFactors = F)
DnDSdata<-na.omit(dfGene2dNdS)

Capacitor<-read.csv('./Fig_3_Results/Capacitor.csv',stringsAsFactors = F)
Potentiator<-read.csv('./Fig_3_Results/Potentiator.csv',stringsAsFactors = F)

Cap_DNDS<-DnDSdata%>%
  dplyr::filter(.,Gene %in% toupper(Capacitor$Gene))%>%
  dplyr::mutate(type="Capacitor")
Pot_DNDS<-DnDSdata%>%
  dplyr::filter(.,Gene %in% toupper(Potentiator$Gene))%>%
  dplyr::mutate(type="Potentiator")

Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)%>%toupper()
Other_DNDS<-DnDSdata%>%
  dplyr::filter(.,!(Gene %in% Co_Gene))%>%
  dplyr::mutate(type="Others")

wilcox.test(Cap_DNDS$dNdS,Other_DNDS$dNdS,alternative = "less")
wilcox.test(Pot_DNDS$dNdS,Other_DNDS$dNdS,alternative = "less")

##draw Fig3_A
Fig3_A<-Cap_DNDS%>%
  rbind(.,Pot_DNDS,Other_DNDS)%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(Mean_DN=mean(dNdS))%>%
  ggplot(aes(x=type,y=Mean_DN,fill=type))+geom_bar(stat="identity")+
  xlab("Gene Type")+ylab("Mean dN/dS of gene")+
  scale_fill_manual(name="Gene Type",
                    labels=c("Capacitor","Others","Potentiator"),
                    values=c("#F8766D","#619CFF","#00BFC4"))+
  scale_x_discrete(limits=c("Capacitor","Potentiator","Others"),
                   labels=c("Stabilizer","Diversifier","Others"))+
  geom_signif(y_position=c(0.078,0.085), xmin=c(2,1), xmax=c(3,3),
              annotation=c("*","0.1"), tip_length=0.1, size=0.8, textsize = 7,vjust = 0.1)+
  theme(legend.position="none")


## Fitness
#setwd("~/phenotypic_heterogeneity/")
Fit_data<-read.table("./Fig_3_Results/fitness_yeast.txt",header = T,stringsAsFactors = F)
Fitness<-Fit_data[,1:2]
Cap_Fit<-Fitness%>%
  dplyr::filter(.,ORF %in% toupper(Capacitor$Gene))%>%
  dplyr::mutate(type="Capacitor")
Pot_Fit<-Fitness%>%
  dplyr::filter(.,ORF %in% toupper(Potentiator$Gene))%>%
  dplyr::mutate(type="Potentiator")
Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)%>%toupper()
Other_Fit<-Fitness%>%
  dplyr::filter(.,!(ORF %in% Co_Gene))%>%
  dplyr::mutate(type="Others")

wilcox.test(Cap_Fit$YPD_TC1,Other_Fit$YPD_TC1,alternative = "less")
wilcox.test(Pot_Fit$YPD_TC1,Other_Fit$YPD_TC1,alternative = "less")
wilcox.test(Cap_Fit$YPD_TC1,Pot_Fit$YPD_TC1,alternative = "less")

##draw Fig3_B
Fig3_B<-Cap_Fit%>%
  rbind(.,Pot_Fit,Other_Fit)%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(Mean_FI=mean(YPD_TC1))%>%
  ggplot(aes(x=type,y=Mean_FI,fill=type))+geom_bar(stat="identity")+
  xlab("Gene Type")+ylab("Fitness upon deletion of gene")+
 scale_fill_manual(name="Gene Type",
                      labels=c("Capacitor","Others","Potentiator"),
                      values=c("#F8766D","#619CFF","#00BFC4"))+
  scale_x_discrete(limits=c("Capacitor","Potentiator","Others"),
                   labels=c("Stabilizer","Diversifier","Others"))+
  geom_signif(y_position=c(1,1.08,0.95), xmin=c(2,1,1), xmax=c(3,3,2),
              annotation=c("***","***","***"), tip_length=0.1, size=0.8, textsize = 7,vjust = 0.1)+
  theme(legend.position="none")


##Fig3_CD

##Trait 70 importance (TI)
load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")
Trait_70_TI<-dfTrait2TI%>%
  dplyr::filter(.,Trait %in% unique(df70trait2Gene_NorDM$Trait))

# Number of Capacitor and Potentiator regulating a trait in 70 Traits
Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)

dfTrait2CapPot_TI<-dfGene2Trait%>%
  dplyr::filter(.,Trait %in%  unique(df70trait2Gene_NorDM$Trait))%>%
  dplyr::filter(.,Gene %in% Co_Gene)%>%
  dplyr::group_by(Trait) %>%
  dplyr::summarise(numAsCap = length(which(Type == "capacitor")),
                   numAsPot = length(which(Type == "potentiator"))) %>%
  merge(.,Trait_70_TI,by="Trait")%>%
  dplyr::arrange(.,TI)

cor.test(dfTrait2CapPot_TI$TI,dfTrait2CapPot_TI$numAsCap,method = "s")
cor.test(dfTrait2CapPot_TI$TI,dfTrait2CapPot_TI$numAsPot,method = "s")

##draw Fig3_C
DATA_Fig3_C<-dfTrait2CapPot_TI%>%
  dplyr::select(.,1,2,4)%>%
  mutate(id=0:69)%>%
  mutate(N_ID=(id%/%7)+1)%>%
  dplyr::group_by(N_ID)%>%
  dplyr::summarise(TI_Mean=mean(TI),Y_Mean=mean(numAsCap),
                   Y_se=sd(numAsCap)/sqrt(length(TI)),
                   Y_Max=Y_Mean+Y_se,Y_Min=Y_Mean-Y_se)

Fig3_C<-DATA_Fig3_C%>%
  ggplot(aes(x=TI_Mean,y=Y_Mean))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=.04)+
  geom_point(color="#F8766D",size=2)+
  xlab(" Mean trait importance ")+
  ylab("Mean number of stabilizer")+
  scale_x_continuous(breaks = seq(0,0.9,0.1))+
  annotate("text",x=.15,y=90,parse = TRUE,label= "'ρ = 0.35,' ~ italic(P) < ~10^-3")

##draw Fig3_D
DATA_Fig3_D<-dfTrait2CapPot_TI%>%
  dplyr::select(.,1,3,4)%>%
  dplyr::mutate(id=0:69)%>%
  dplyr::mutate(N_ID=(id%/%7)+1)%>%
  dplyr::group_by(N_ID)%>%
  dplyr::summarise(TI_Mean=mean(TI),Y_Mean=mean(numAsPot),
                   Y_se=sd(numAsPot)/sqrt(length(TI)),
                   Y_Max=Y_Mean+Y_se,Y_Min=Y_Mean-Y_se)

Fig3_D<- DATA_Fig3_D%>%
  ggplot(aes(x=TI_Mean,y=Y_Mean))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=.04)+
  geom_point(color="#00BFC4",size=2)+
  xlab(" Mean trait importance ")+
  ylab("Mean number of diversifier")+
  scale_x_continuous(breaks = seq(0,0.9,0.1))+
  annotate("text",x=.12,y=90,parse = TRUE,label= "'ρ = -0.3,' ~ italic(P) < 0.012")


##Fig3_EF

dfGene2dNdS<-read.csv("./Fig_3_Results/dNdS.csv",stringsAsFactors = F)
DnDSdata<-na.omit(dfGene2dNdS)
Traits70<-unique(df70trait2Gene_NorDM$Trait)
Co_Gene<-union(Capacitor$Gene,Potentiator$Gene)

df70Trait_Gene<-dfGene2Trait%>%
  dplyr::filter(.,Trait %in% Traits70)%>%
  dplyr::filter(.,Gene %in% Co_Gene)%>%
  dplyr::mutate(AGene=toupper(Gene))%>%
  merge(DnDSdata,.,by.x="Gene",by.y="AGene")%>%
  merge(dfTrait2TI,by="Trait")%>%
  dplyr::select(.,1,2,3,5,6)

CapDATA<-filter(df70Trait_Gene,Type=="capacitor" )  
PotDATA<-filter(df70Trait_Gene,Type=="potentiator") 

Test_Cap_DNDS<-mclapply(1:70,function(k){
  A_Trait<-filter(CapDATA,Trait %in% Traits70[k])
  A_Trait_Own_Cap<-A_Trait[,c(2,3)] #chose overlap capacitor and DNDS
  if(length( A_Trait_Own_Cap$Gene)<2 |length(A_Trait_Own_Cap$dNdS)==0)
  {new<-data.frame(DATA="NA")}else{
    new<-data.frame(Trait=Traits70[k],Mean_DNDS=mean(A_Trait_Own_Cap$dNdS),
                    TI=A_Trait$TI[1])
  }
  
  return(new)
},mc.cores = 10)%>%rbind.fill()
Test_Pot_DNDS<-mclapply(1:70,function(k){
  A_Trait<-filter(PotDATA,Trait %in% Traits70[k])
  A_Trait_Own_Pot<-A_Trait[,c(2,3)] #chose overlap potentiator and DNDS
  
  if(length(A_Trait_Own_Pot$dNdS)==0)
  {new_1<-data.frame(DATA="NA")}else{
    new_1<-data.frame(Trait=Traits70[k],Mean_DNDS=mean(A_Trait_Own_Pot$dNdS),
                      TI=A_Trait$TI[1])
  }
  
  return(new_1)
  
},mc.cores = 10)%>%rbind.fill()

Test_Cap_DNDS$Trait<-as.character(Test_Cap_DNDS$Trait)
Test_Pot_DNDS$Trait<-as.character(Test_Pot_DNDS$Trait)

Test_Cap_DNDS<-Test_Cap_DNDS%>%
  dplyr::filter(.,nchar(Trait)!=0)%>%
  dplyr::filter(., !(Trait %in% c("A103-C","A123-A","C13-A"))) #Remove outliers

cor.test(Test_Cap_DNDS$Mean_DNDS,Test_Cap_DNDS$TI,method = "s")
cor.test(Test_Pot_DNDS$Mean_DNDS,Test_Pot_DNDS$TI,method = "s")

##draw Fig3_E

Fig3_E_DATA<-Test_Cap_DNDS%>%
  dplyr::arrange(.,TI)%>%
  dplyr::mutate(id=0:(nrow(Test_Cap_DNDS)-1),
                Group_ID=(id%/%7)+1)%>%
  dplyr::group_by(Group_ID)%>%
  dplyr::summarise(TI_Mean=mean(TI),Y_Mean=mean(Mean_DNDS),
                   Y_se=sd(Mean_DNDS)/sqrt(length(TI)),
                   Y_Max=Y_Mean+Y_se,Y_Min=Y_Mean-Y_se)

Fig_3_E<- Fig3_E_DATA%>%
  ggplot(aes(x=TI_Mean,y=Y_Mean))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=.04)+
  geom_point(color="#F8766D",size=2)+
  xlab(" Mean trait importance ")+
  ylab("Mean dN/dS of stabilizer")+
  scale_x_continuous(breaks = seq(0,0.9,0.1))+
  annotate("text",x=.15,y=0.08,parse = TRUE,label= "'ρ = 0.056,' ~ italic(P) == 0.6642")


##draw Fig3_F
Fig3_F_DATA<-Test_Pot_DNDS%>%
  dplyr::arrange(.,TI)%>%
  dplyr::mutate(id=0:(nrow(Test_Pot_DNDS)-1),
                Group_ID=(id%/%7)+1)%>%
  dplyr::group_by(Group_ID)%>%
  dplyr::summarise(TI_Mean=mean(TI),Y_Mean=mean(Mean_DNDS),
                   Y_se=sd(Mean_DNDS)/sqrt(length(TI)),
                   Y_Max=Y_Mean+Y_se,Y_Min=Y_Mean-Y_se)

Fig3_F<-Fig3_F_DATA%>%
  ggplot(aes(x=TI_Mean,y=Y_Mean))+
  geom_errorbar(aes(ymin=Y_Min,ymax=Y_Max),width=.04)+
  geom_point(color="#00BFC4",size=2)+
  xlab(" Mean trait importance ")+
  ylab("Mean dN/dS of diversifier")+
  scale_x_continuous(breaks = seq(0,0.9,0.1))+
  annotate("text",x=.15,y=0.08,parse = TRUE,label="'ρ = -0.28,' ~ italic(P) < 0.03")



##Draw  DATA save in Figure_3.Rdata


