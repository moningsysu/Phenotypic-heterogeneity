library(readr);
library(dplyr);
library(parallel);
library(ggsignif);

##Draw  DATA save in Figure_3.Rdata

load("~/phenotypic_heterogeneity/dNdS.RData")
load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")

Traits70<-unique(df70trait2Gene_NorDM$Trait)
Co_Gene <- union(Capacitor$Gene,Potentiator$Gene)%>%toupper()
Gene4718<- toupper(unique(df220Trait_4718Gene_DM$Gene))

##
## 1.  average dN/dS of each gene in the five comparison
##
df.mean.dNdS <- df.dNdS %>%
  dplyr::filter(ORF %in% toupper(Gene4718)) %>%
  dplyr::group_by(ORF) %>%
  dplyr::summarise(mean_dN_dS = mean(dN_dS,na.rm=T));

vCap <- df.mean.dNdS %>%
  filter(ORF %in% toupper(Capacitor$Gene)) %>%
  getElement("mean_dN_dS")

vPot <- df.mean.dNdS %>%
  filter(ORF %in% toupper(Potentiator$Gene)) %>%
  getElement("mean_dN_dS")

vOther <- df.mean.dNdS %>%
  filter(ORF %in% toupper(Gene4718)) %>%
  filter(!(ORF %in% toupper(c(Capacitor$Gene,Potentiator$Gene)) )) %>%
  getElement("mean_dN_dS")

wilcox.test(vCap,vOther,alternative = "less");
wilcox.test(vPot,vOther,alternative = "less");
wilcox.test(vPot,vCap,alternative = "less");

## plot Fig. 3A using vCap, vPot, vOther
##
toPlot.3A <- data.frame(dnds = c(vCap,vPot,vOther),
                        type = ordered(
                          rep(c("Stabilizer","Diversifier","Other"),
                              c(length(vCap),length(vPot),length(vOther)) ),
                          levels=c("Stabilizer","Diversifier","Other")))

##draw Fig3_A
Fig3_A <-  toPlot.3A %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(mean = mean(dnds),
                   sd = sd(dnds),
                   se = sd(dnds) / sqrt(length(dnds))) %>%
  ggplot(aes(x=type,y=mean,fill=type))+
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2) +
  scale_fill_manual(name="Gene type",
                    labels=c("Stabilizer","Diversifier","Other"),
                    values=c("#F8766D","#00BFC4","#619CFF")) +
  scale_color_manual(name="Gene Type",
                     labels=c("Stabilizer","Other","Diversifier"),
                     values=c("#F8766D","#619CFF","#00BFC4")) +
  scale_x_discrete("Gene type") +
  scale_y_continuous("Mean dN/dS of genes",limits=c(0,0.12)) +
  geom_signif(color="black",
              y_position=c(0.105,0.115), xmin=c(2,1), xmax=c(3,3),
              annotation=c("**","*"), tip_length=0.1, size=0.8, 
              textsize = 7,vjust = 0.5)+
  theme(legend.position = "none")

##
## Fig3B
##
Fit_data<-read.table("~/phenotypic_heterogeneity/Fig_3_Results/fitness_yeast.txt",header = T,stringsAsFactors = F)
fCap <- Fit_data %>%
  filter(ORF %in% toupper(Capacitor$Gene)) %>%
  getElement("YPD_TC1")
fPot <- Fit_data %>%
  filter(ORF %in% toupper(Potentiator$Gene)) %>%
  getElement("YPD_TC1")
fOther <- Fit_data %>%
  filter(ORF %in% toupper(Gene4718)) %>%
  filter(!(ORF %in% toupper(c(Capacitor$Gene,Potentiator$Gene)) )) %>%
  getElement("YPD_TC1")
wilcox.test(fCap,fOther,alternative = "less");
wilcox.test(fPot,fOther,alternative = "less");
wilcox.test(fCap,fPot,alternative = "less");
toPlot.3b <- data.frame(fit = c(fCap,fPot,fOther),
                        type = ordered(
                          rep(c("Stabilizer","Diversifier","Other"),
                              c(length(fCap),length(fPot),length(fOther)) ),
                          levels=c("Stabilizer","Diversifier","Other")))

Fig3_B<-toPlot.3b%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(Mean_FI=mean(fit),
                   sd_FI=sd(fit),
                   se_FI=sd(fit)/sqrt(length(fit)))%>%
  ggplot(aes(x=type,y=Mean_FI,fill=type))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=Mean_FI-se_FI,ymax=Mean_FI+se_FI),width=0.2) +
  xlab("Gene type")+ylab("Fitness upon deletion of gene  ")+
  scale_fill_manual(name="Gene Type",
                    labels=c("Stabilizer","Diversifier","Other"),
                    values=c("#F8766D","#00BFC4","#619CFF"))+
  scale_x_discrete(limits=c("Stabilizer","Diversifier","Other"),
                   labels=c("Stabilizer","Diversifier","Other"))+
  geom_signif(y_position=c(1,1.08), xmin=c(2,1), xmax=c(3,3),
              annotation=c("***","***"), tip_length=0.1, size=0.8, textsize = 7,vjust = 0.5)+
  theme(legend.position="none")

##
## Fig3C/D

#load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")
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

Fig3_C<-dfTrait2CapPot_TI%>%
  dplyr::select(.,1,2,4)%>%
  ggplot(aes(x=TI,y=numAsCap))+
  geom_point(color="#F8766D",size=1.5)+
  geom_smooth(method="lm",se = F,color="grey")+
  xlab("Trait importance    ")+
  ylab("Number of stabilizer")+
  scale_x_log10()+
  annotate("text",x=0.05,y=120,parse = TRUE,label= "rho ~ ' = 0.35,'")+
  annotate("text",x=0.05,y=105,parse = TRUE,label= "~ italic(P) < ~10^-3")

## draw Fig3_D
## compare top/bottom 20 traits
tryProb <- 0.3; ## so that the top/bottom 20 traits were compared
nGene.cap.lowTi <- dfTrait2CapPot_TI %>%
  filter(TI < quantile(TI,p=tryProb)) %>%
  getElement("numAsCap")
nGene.cap.highTi <- dfTrait2CapPot_TI %>%
  filter(TI > quantile(TI,p=1-tryProb)) %>%
  getElement("numAsCap")
wilcox.test(nGene.cap.lowTi,nGene.cap.highTi); # P = 0.02
toPlot.3d <- data.frame(nGene = c(nGene.cap.lowTi,nGene.cap.highTi),
                        type = ordered(
                          rep(c("low","high"),
                              c(length(nGene.cap.lowTi),length(nGene.cap.highTi)) ),
                          
                          levels=c("low","high")));

Fig3_D<-toPlot.3d%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(meanN=mean(nGene),sdN=sd(nGene))%>%
  ggplot(aes(x=type,y=meanN))+
  geom_bar(stat="identity",fill="#F8766D")+
  geom_errorbar(aes(ymax=meanN+sdN/sqrt(20),ymin=meanN-sdN/sqrt(20)),
                color="black",width=0.2)+
  xlab("")+
  scale_y_continuous("\nMean number of stabilizer   ",limits = c(0,73))+
  scale_x_discrete(limits=c("low","high"),
                   labels=c("The 20 least\nimportant traits",
                            "The 20 most\nimportant traits"))+
  geom_signif(y_position=c(68), xmin=c(1), xmax=c(2),
              annotation= "*",
              parse=F,
              tip_length=0.1, size=0.8, textsize = 7)+
  theme(legend.position="none",axis.text.x=element_text(angle=40,vjust=1,hjust=1))


## Fig3 E/F

## draw Fig3 E
Fig3_E<-dfTrait2CapPot_TI%>%
  dplyr::select(.,1,3,4)%>%
  ggplot(aes(x=TI,y=numAsPot))+
  geom_point(color="#00BFC4",size=1.5)+
  geom_smooth(method="lm",se = F,color="grey")+
  xlab("Trait importance    ")+
  ylab("Number of diversifier")+
  scale_x_log10()+
  annotate("text",x=0.04,y=145,parse = TRUE,label= "rho ~ ' = -0.30,'")+
  annotate("text",x=0.04,y=130,parse = TRUE,label= "~ italic(P) < 0.02")

## draw Fig3 F
nGene.pot.lowTi <- dfTrait2CapPot_TI %>%
  filter(TI < quantile(TI,p=tryProb)) %>%
  getElement("numAsPot")
nGene.pot.highTi <- dfTrait2CapPot_TI %>%
  filter(TI > quantile(TI,p=1-tryProb)) %>%
  getElement("numAsPot")
wilcox.test(nGene.pot.lowTi,nGene.pot.highTi); ## P = 0.0031
toPlot.3f <- data.frame(nGene = c(nGene.pot.lowTi,nGene.pot.highTi),
                        type = ordered(
                          rep(c("low","high"),
                              c(length(nGene.pot.lowTi),length(nGene.pot.highTi)) ),
                          levels=c("low","high")));

Fig3_F<-toPlot.3f%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(meanN=mean(nGene),sdN=sd(nGene))%>%
  ggplot(aes(x=type,y=meanN))+
  geom_bar(stat="identity",fill="#00BFC4")+
  geom_errorbar(aes(ymax=meanN+sdN/sqrt(20),ymin=meanN-sdN/sqrt(20)),
                color="black",width=0.2)+
  scale_y_continuous("\nMean number of diversifier   ",limits = c(0,77))+
  scale_x_discrete("",
                   limits=c("low","high"),
                   labels=c("The 20 least\nimportant traits",
                            "The 20 most\nimportant traits"))+
  geom_signif(y_position=c(71), xmin=c(1), xmax=c(2),
              annotation= "**",
              parse=F,
              tip_length=0.1, size=0.8, textsize = 7)+
  theme(legend.position="none",axis.text.x=element_text(angle=40,vjust=1,hjust=1))




## Fig3 G/I

##load("~/phenotypic_heterogeneity/Original_Data/Original.Rdata")
Co_Gene <- union(Capacitor$Gene,Potentiator$Gene)%>%toupper()
dfGene2Trait <- dfGene2Trait %>% mutate(Gene = toupper(Gene));
df70Trait_Gene <- dfGene2Trait%>%
  dplyr::filter(Trait %in% Traits70)%>%
  dplyr::filter(Gene %in% Co_Gene)%>%
  mutate(Gene = toupper(Gene)) %>%
  merge(df.mean.dNdS,.,by.x="ORF",by.y="Gene")%>%
  merge(dfTrait2TI,by="Trait");

CapDATA<-filter(df70Trait_Gene,Type=="capacitor" )  ;
PotDATA<-filter(df70Trait_Gene,Type=="potentiator") ;

Cap_DNDS<-mclapply(1:70,function(k){
  A_Trait<-filter(CapDATA,Trait %in% Traits70[k])
  A_Trait_Own_Cap<-A_Trait[,c("ORF","mean_dN_dS")] #chose overlap capacitor and DNDS
  if(length( A_Trait_Own_Cap$ORF)<2 |length(A_Trait_Own_Cap$`mean_dN_dS`)==0){
    new<-data.frame(NULL)
  }
  else    {
    new<-data.frame(Trait=Traits70[k],Mean_DNDS=mean(A_Trait_Own_Cap$`mean_dN_dS`),
                    TI=A_Trait$TI[1])
  }
  
  return(new)
},mc.cores = 10)%>%rbind.fill()

Pot_DNDS<-mclapply(1:70,function(k){
  A_Trait<-filter(PotDATA,Trait %in% Traits70[k])
  A_Trait_Own_Pot<-A_Trait[,c("ORF","mean_dN_dS")] #chose overlap potentiator and DNDS
  
  if(length(A_Trait_Own_Pot$ORF)==0)
  {new_1<-data.frame(NULL)}else{
    new_1<-data.frame(Trait=Traits70[k],Mean_DNDS=mean(A_Trait_Own_Pot$`mean_dN_dS`),
                      TI=A_Trait$TI[1])
  }
  
  return(new_1)
  
},mc.cores = 10)%>%rbind.fill()

Cap_DNDS$Trait<-as.character(Cap_DNDS$Trait)
Pot_DNDS$Trait<-as.character(Pot_DNDS$Trait)

cor.test(Cap_DNDS$Mean_DNDS,Cap_DNDS$TI,method = "s")
cor.test(Pot_DNDS$Mean_DNDS,Pot_DNDS$TI,method = "s")

##
## plot Fig. 3G/H/I/J based on the correlation above

##draw Fig3 G
Fig3_G <- Cap_DNDS %>%
  ggplot(aes(x=TI,y=Mean_DNDS))+
  geom_point(color="#F8766D",size=1.5)+
  geom_smooth(method="lm",se = F,color="grey")+
  xlab("Trait importance     ")+
  ylab("dN/dS of stabilizer")+
  scale_x_log10()+
  annotate("text",x=0.04,y=0.1,parse = TRUE,label= "rho ~ ' = -0.11,'")+
  annotate("text",x=0.04,y=0.09,parse = TRUE,label= "~ italic(P) == 0.39")



## draw Fig3 H
## compare top/bottom 20 traits
tryProb <- 0.3; ## so that the top/bottom 20 traits were compared
dnds.cap.lowTi <- Cap_DNDS %>%
  filter(TI < quantile(TI,p=tryProb)) %>%
  getElement("Mean_DNDS")
dnds.cap.highTi <- Cap_DNDS %>%
  filter(TI > quantile(TI,p=1-tryProb)) %>%
  getElement("Mean_DNDS")
wilcox.test(dnds.cap.lowTi,dnds.cap.highTi); ## P = 0.43
toPlot.3h <- data.frame(dnds = c(dnds.cap.lowTi,dnds.cap.highTi),
                        type = ordered(
                          rep(c("low","high"),
                              c(length(dnds.cap.lowTi),length(dnds.cap.highTi)) ),
                          levels=c("low","high")))
Fig3_H<-toPlot.3h%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(meanDnds=mean(dnds),sdDnds=sd(dnds))%>%
  ggplot(aes(x=type,y=meanDnds))+
  geom_bar(stat="identity",fill="#F8766D")+
  geom_errorbar(aes(ymax=meanDnds+sdDnds/sqrt(20),ymin=meanDnds-sdDnds/sqrt(20)),
                color="black",width=0.2)+
  xlab("")+
  scale_y_continuous("\nMean dN/dS of stabilizer    ",limits = c(0,0.11))+
  scale_x_discrete(limits=c("low","high"),
                   labels=c("The 20 least\nimportant traits",
                            "The 20 most\nimportant traits"))+
  geom_signif(y_position=c(0.1), xmin=c(1), xmax=c(2),
              annotation= "italic(P) ~ \" = 0.43\"",
              parse=T,
              tip_length=0.1, size=0.8)+
  theme(legend.position="none",axis.text.x=element_text(angle=40,vjust=1,hjust=1))


##draw Fig3 I


Fig3_I <- Pot_DNDS %>%
  ggplot(aes(x=TI,y=Mean_DNDS))+
  geom_point(color="#00BFC4",size=1.5)+
  geom_smooth(method="lm",se = F,color="grey")+
  xlab("trait importance     ")+
  ylab(" dN/dS of diversifier")+
  scale_x_log10()+
  scale_y_continuous(limits=c(0.04,0.09))+
  annotate("text",x=0.05,y=0.085,parse = TRUE,label="rho ~' = -0.26,'")+
  annotate("text",x=0.05,y=0.080,parse = TRUE,label="~ italic(P) < 0.04")

## draw Fig3 J
## compare top/bottom 20 traits
dnds.pot.lowTi <- Pot_DNDS %>%
  filter(TI < quantile(TI,p=tryProb)) %>%
  getElement("Mean_DNDS")
dnds.pot.highTi <- Pot_DNDS %>%
  filter(TI > quantile(TI,p=1-tryProb)) %>%
  getElement("Mean_DNDS")
wilcox.test(dnds.pot.lowTi,dnds.pot.highTi)## P = 0.0067
toPlot.3j <- data.frame(dnds = c(dnds.pot.lowTi,dnds.pot.highTi),
                        type = ordered(
                          rep(c("low","high"),
                              c(length(dnds.cap.lowTi),length(dnds.cap.highTi)) ),
                          levels=c("low","high")))
Fig3_J <-toPlot.3j%>%
  dplyr::group_by(type)%>%
  dplyr::summarise(meanDnds=mean(dnds),sdDnds=sd(dnds))%>%
  ggplot(aes(x=type,y=meanDnds))+
  geom_bar(stat="identity",fill="#00BFC4")+
  geom_errorbar(aes(ymax=meanDnds+sdDnds/sqrt(20),ymin=meanDnds-sdDnds/sqrt(20)),
                color="black",width=0.2)+
  xlab("")+
  scale_y_continuous("\nMean dN/dS of diversifier    ",limits = c(0,0.09))+
  scale_x_discrete(limits=c("low","high"),
                   labels=c("The 20 least\nimportant traits",
                            "The 20 most\nimportant traits"))+
  geom_signif(y_position=c(0.08), xmin=c(1), xmax=c(2),
              annotation= "**",
              parse=F,
              tip_length=0.1, size=0.8, textsize = 7)+
  theme(legend.position="none",axis.text.x=element_text(angle=40,vjust=1,hjust=1))



##Draw  DATA save in Figure_3.Rdata
save(df.dNdS,Fit_data,dfTrait2CapPot_TI,Cap_DNDS,Pot_DNDS,
     file = "~/12_5/Figure_3.Rdata")



