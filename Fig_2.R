library(ggplot2);
library(openxlsx);
library(reshape2);
library(plyr);
library(dplyr);
library(scales);
library(cowplot);

setwd("~/phenotypic_heterogeneity/Fig_2_Results/")
dfRawdm <- read.xlsx("./TabS1.rawDMs.xlsx",sheet=1,startRow=3);
dfRoles <- read.xlsx("./TabS2.roles.xlsx",sheet=1,startRow=2);
dfLinks <- read.xlsx("./TabS3.network.xlsx",sheet=1,startRow=2);


## add NorDm to the links of the network
dfRawdm[,-1] <- dfRawdm[,-1] %>% apply(2,scale);
dfRawdm.melted <- melt(dfRawdm,id.vars="Gene");
colnames(dfRawdm.melted) <- c("Gene","Trait","standardDM");
dfLinks.withDm <- dfLinks %>% merge(dfRawdm.melted,by=c("Gene","Trait"))



## Trait  view
## Number of genes regulating a trait

dfLinks.eachTrait <- dfLinks.withDm %>%
  merge(dfRoles,by.x="Gene",by.y="Gene") %>%
  filter( ((Type_of_regulation == "Potentiation") & (Type == "Potentiator")) |
            ((Type_of_regulation == "Capacitance") & (Type == "Capacitor"))  ) %>%
  group_by(Trait,Type_of_regulation) %>%
  dplyr::summarise(numReg = length(Trait),avgAbsDm = mean(abs(standardDM)));

nReg.eachTrait <- split(dfLinks.eachTrait$numReg,dfLinks.eachTrait$Type_of_regulation);
lapply(nReg.eachTrait,median);
wilcox.test(nReg.eachTrait$Capacitance,nReg.eachTrait$Potentiation);

##Fig2_A
dfLinks.eachTrait$numReg %>% quantile(p=seq(0,1,0.05)); ## decide the breaks for the plot
toPlot.dfLinks.eachTrait <- dfLinks.eachTrait %>%
  mutate(numRegLvl = cut(numReg,breaks=c(seq(0,120,10),Inf),ordered_result=T)) %>%
  group_by(numRegLvl,Type_of_regulation) %>%
  dplyr::summarise(cnt = length(Trait)) %>%
  ungroup() %>%
  group_by(Type_of_regulation) %>%
  dplyr::mutate(frac = cnt/sum(cnt));

plot.trait.nReg <- toPlot.dfLinks.eachTrait %>%
  ggplot(aes(x=numRegLvl,fill=Type_of_regulation,y=frac)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5,preserve="single")) +
  scale_x_discrete("No. of genes regulating a trait") +
  scale_y_continuous("Fraction of traits") +
  scale_fill_discrete("Gene Type ",labels=c("Capacitance"="Stabilizer","Potentiation"="Diversifier")) +
  theme_classic() +
  theme(legend.position = c(1,1),legend.justification = c(1,1),
        axis.text.x = element_text(angle = 45, hjust = 1));

## Fig2_B  Average effect size over a trait
avgAbsDm.eachTrait <- split(dfLinks.eachTrait$avgAbsDm,dfLinks.eachTrait$Type_of_regulation);
lapply(avgAbsDm.eachTrait,median);
wilcox.test(avgAbsDm.eachTrait$Capacitance,avgAbsDm.eachTrait$Potentiation);

dfLinks.eachGene$avgAbsDm %>% quantile(p=seq(0,1,0.05)); ## decide the breaks for the plot
toPlot.effSize.dfLinks.eachTrait <- dfLinks.eachTrait %>%
  mutate(avgAbsDmLvl = cut(avgAbsDm,breaks=c(seq(0,6,0.5),Inf),ordered_result=T)) %>%
  group_by(avgAbsDmLvl,Type_of_regulation) %>%
  dplyr::summarise(cnt = length(Trait)) %>%
  ungroup() %>%
  group_by(Type_of_regulation) %>%
  dplyr::mutate(frac = cnt/sum(cnt));

plot.trait.effSize <- toPlot.effSize.dfLinks.eachTrait %>%
  ggplot(aes(x=avgAbsDmLvl,fill=`Type_of_regulation`,y=frac)) +
  geom_bar(stat="identity",position=position_dodge(preserve="single",width=0.5)) +
  scale_x_discrete("Average effect size over a trait (absolute standardized DM)") +
  scale_y_continuous(("Fraction of traits")) +
  scale_fill_discrete("Gene Type ",labels=c("Capacitance"="Stabilizer","Potentiation"="Diversifier")) +
  theme_classic() +
  theme(legend.position=c(1,1),legend.justification=c(1,1),
        axis.text.x = element_text(angle = 45, hjust = 1));




##Fig2_C
##  Number of traits regulated
dfLinks.eachGene <- dfLinks.withDm %>%
  merge(dfRoles,by.x="Gene",by.y="Gene") %>%
  filter( ((Type_of_regulation == "Potentiation") & (Type == "Potentiator")) |
            ((Type_of_regulation == "Capacitance") & (Type == "Capacitor"))  ) %>%
  dplyr::group_by(Gene,Type_of_regulation,Type) %>%
  dplyr::summarise(numReg = length(Gene),avgAbsDm = mean(abs(standardDM)));

nReg.eachGene <- split(dfLinks.eachGene$numReg,dfLinks.eachGene$Type_of_regulation);
lapply(nReg.eachGene,median);
wilcox.test(nReg.eachGene$Capacitance,nReg.eachGene$Potentiation);


dfLinks.eachGene$numReg %>% quantile(p=seq(0,1,0.05)); ## decide the breaks for the plot
toPlot.dfLinks.eachGene <- dfLinks.eachGene %>%
  mutate(numRegLvl = cut(numReg,breaks=c(seq(0,14,1),Inf),ordered_result=T)) %>%
  group_by(numRegLvl,Type_of_regulation) %>%
  dplyr::summarise(cnt = length(Gene)) %>%
  ungroup() %>%
  group_by(Type_of_regulation) %>%
  dplyr::mutate(frac = cnt/sum(cnt));

plot.gene.nReg <- toPlot.dfLinks.eachGene %>%
  ggplot(aes(x=as.numeric(numRegLvl),fill=Type_of_regulation,y=frac)) +
  geom_bar(stat="identity",position=position_dodge(width=0.5)) +
  scale_x_continuous("No. of traits regulated by a gene",
                     breaks=c(1:15),
                     labels=expression(1,2,3,4,5,6,7,8,9,10,11,12,13,14,"" > 14)) +
  scale_y_continuous("Fraction of genes") +
  scale_fill_discrete("Gene Type ",labels=c("Capacitance"="Stabilizer","Potentiation"="Diversifier")) +
  theme_classic()+
  theme(legend.position = c(1,1),legend.justification = c(1,1));



## Fig2_D Average effect size of a gene
avgAbsDm.eachGene <- split(dfLinks.eachGene$avgAbsDm,dfLinks.eachGene$Type_of_regulation);
lapply(avgAbsDm.eachGene,median);
wilcox.test(avgAbsDm.eachGene$Capacitance,avgAbsDm.eachGene$Potentiation);

dfLinks.eachGene$avgAbsDm %>% quantile(p=seq(0,1,0.05)); ## decide the breaks for the plot
toPlot.effSize.dfLinks.eachGene <- dfLinks.eachGene %>%
  mutate(avgAbsDmLvl = cut(avgAbsDm,breaks=c(seq(0,5.2,0.4),Inf),ordered_result=T)) %>%
  group_by(avgAbsDmLvl,Type_of_regulation) %>%
  dplyr::summarise(cnt = length(Gene)) %>%
  ungroup() %>%
  group_by(Type_of_regulation) %>%
  dplyr::mutate(frac = cnt/sum(cnt));

plot.gene.effSize <- toPlot.effSize.dfLinks.eachGene %>%
  ggplot(aes(x=avgAbsDmLvl,fill=`Type_of_regulation`,y=frac)) +
  geom_bar(stat="identity",position=position_dodge(preserve="single",width=0.5)) +
  scale_x_discrete("Average effect size of a gene (absolute standardized DM)") +
  scale_y_continuous(("Fraction of genes")) +
  scale_fill_discrete("Gene Type ",labels=c("Capacitance"="Stabilizer","Potentiation"="Diversifier")) +
  theme_classic() +
  theme(legend.position=c(1,1),legend.justification=c(1,1),
        axis.text.x = element_text(angle = 45, hjust = 1))

##Draw  DATA save in Figure_2.Rdata


