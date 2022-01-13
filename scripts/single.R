#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0)
  stop("Provide the harmony.txt files\n", call.=FALSE)

library(vroom)
library(dplyr)
library(ggplot2)

data = NULL

for (p in args) {
  data_a = vroom(p,delim = ' ')
  data_a$dataset = p
  data = bind_rows(data,data_a)
}

qv_summary = data %>% group_by(ec,dataset) %>% filter(ec < 30)  %>% summarise(numBases = sum(alnlen), numMM = sum(mismatch), numM = sum(match), numIE = sum(ins_events), numDE = sum(del_events), qv = -10*log10(1-(1+numM)/(1+numM+numMM+numIE+numDE)))
# qv_summary = data %>% group_by(ec,dataset) %>% filter(ec < 30) %>% summarise(numBases = sum(alnlen), TotalBases = sum(match + mismatch + ins), error = sum(mismatch + del + ins), qv = ifelse(error==0,ceiling(-10*log10(1/TotalBases)), -10*log10(error/TotalBases)))
g1 = ggplot(qv_summary) +
  stat_smooth(aes(ec,qv,fill=dataset,col=dataset), alpha = 0.1) +
  geom_hline(yintercept = 30, lty = 3, alpha = 0.75) +
  geom_vline(xintercept = 15, lty = 3, alpha = 0.75) +
  geom_point(aes(ec,qv,size=numBases,col=dataset), alpha = 0.75) +
  coord_cartesian(ylim=c(20,40), xlim=c(2,20))+
  ylab("Mean Gap-Compressed Identity Phred")+
  xlab("Effective coverage")+
  theme_minimal()+
  scale_size_area(breaks=c(1 %o% 10^(1:20)))+
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="right",
        axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14), legend.title = element_text(size=14,face="bold"))+
  scale_x_continuous(breaks=c(3,5,10,15,20))

ggsave("harmony.pdf",g1,width=25,height=15,dpi=200,units="cm",limitsize = FALSE)
