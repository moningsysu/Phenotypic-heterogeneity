library(plyr);
library(dplyr);
library(scales);
library(cowplot);
library(ggplot2);
library(reshape2);
#(Fig. 1A) 
#For each strain, one nonessential gene was deleted, and various morphological traits were measured for >100 cells at each of three cell cycle stages (Fig. 1A)

One_Gene<-read.csv("~/phenotypic_heterogeneity/Fig_1_Results/yal012w.csv",stringsAsFactors = F)
Draw<-One_Gene%>%
  dplyr::select(cell_id,Dgroup,C11.2,C104)%>%
  dplyr::filter(Dgroup %in% c("A","A1","B","C"))%>%
  dplyr::mutate(Type=ifelse(Dgroup=="A","A",
                            ifelse(Dgroup=="C","C","A1B")))

#Trait C104：Short axis length in mother
Fig_1_A_1<-Draw%>%
  dplyr::select(cell_id,Type,C104)%>%
  reshape2::melt(id.var=c("cell_id","Type"))%>%
  dplyr::filter(value != -1)%>%
  ggplot(aes(x=cell_id,y=value,shape=Type))+
  geom_point(size=1.5,alpha=0.8,color="#FF9933")+
  xlab("Individual cells")+
  ylab("Short-axis length of mother cell")+
  scale_shape_manual(name="Nuclear stage",values=c(19,4,1))


#Trait C11.2 Size of daughter cell
Fig_1_A_2<-Draw%>%
  dplyr::select(cell_id,Type,C11.2)%>%
  reshape2::melt(id.var=c("cell_id","Type"))%>%
  dplyr::filter(value != -1)%>%
  ggplot(aes(x=cell_id,y=value,shape=Type))+
  geom_point(size=1.5,alpha=0.8,color="#006241")+
  xlab("Individual cells")+
  ylab("Size of daughter cell")+
  scale_shape_manual(name="Nuclear stage",values=c(4,1))



#(Fig.1B,Fig.1C):For each trait and each strain,we calculated the coefficient of variation (CV) and its deviation from the median CV of strains with similar mean trait value (DM) 
One_Trait<-read.csv("~/phenotypic_heterogeneity/Fig_1_Results/C106_A1B.csv",stringsAsFactors = F)
Fig_1_B<-One_Trait%>%
  ggplot(aes(x=Mean,y=Cv))+
  geom_point(size=1.5,shape=1,alpha=0.8)+
  ylim(c(0.2,1))+
  xlab("Mean trait value of each mutant")+
  ylab("CV of trait value of each mutant")

Fig_1_C<-One_Trait%>%
  ggplot(aes(x=Mean,y=DM))+
  geom_point(size=1.5,shape=1,alpha=0.8)+
  xlab("Mean trait value of each mutant")+
  ylab("DM of trait value of each mutant")


#(Fig. 1D):we respectively averaged the bottom and top 20 (of 70) standardized DM for each gene 
Draw_Gene<-read.csv("~/phenotypic_heterogeneity/Fig_1_Results/Standardied DM.csv",stringsAsFactors = F)

Fig_1_D<-Draw_Gene%>%
  dplyr::arrange(NorDM)%>%
  plyr::mutate(id=1:70)%>%
  plyr::mutate(TYPE=ifelse(id<21,"Low",
                           ifelse(id>50,"High","Median")))%>%
  ggplot(aes(x=id,y=NorDM,color=TYPE,shape=TYPE))+
  geom_point(size=1.2)+
  scale_x_continuous(breaks = seq(0,70,10))+
  scale_color_manual(breaks = c("High","Low","Median"),
                     values=c("#F8766D","#00BFC4","black"))+
  scale_shape_manual(values=c(19,19,1))+
  xlab("Individual traits")+
  ylab(" Standardized DM \n upon deletion of a gene")


#(Fig. 1E):Comparing observed DMs of an YKO over a trait with 1,000 mock DMs calculated by random samples of wildtype cells 
load("~/phenotypic_heterogeneity/Fig_1_Results/Random_DM.Rdata")

#get Pvalue_Data
Fig_1_E<-Pvalue_Data%>%
  ggplot(aes(x=Random_DM))+
  geom_density()+
  xlim(c(-0.16,0.4))+
  xlab("DM of random wildtype samples")+
  ylab("Frequency")+
  scale_y_continuous(breaks=seq(0,6,2),labels = c("0","20%","40%","60%"))


###########################################################################

###we generated a bipartite network between 673 genes and 70 morphological traits, 
##with links representing potentiation or capacitance relationship between them 
##(Fig. 1F)


##Draw  DATA save in Figure_1.Rdata

