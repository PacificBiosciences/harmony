#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0)
  stop("Provide the profiles.txt files\n", call.=FALSE)

library(vroom)
library(dplyr)
library(ggplot2)

data = NULL

for (p in args) {
  data_a = vroom(p,delim = ' ')
  data_a$dataset = p
  data_a = data_a %>% mutate(TotalBases = (match + mismatch + ins) +1 , error = (mismatch + del + ins) + 1, mi = -10*log10(error/TotalBases))
  data = bind_rows(data,data_a)
}

qv_summary = data %>% group_by(ec,dataset) %>% filter(ec < 60) %>% summarise(numBases = sum(alnlen), TotalBases = sum(match + mismatch + ins), error = sum(mismatch + del + ins), qv = ifelse(error==0,ceiling(-10*log10(1/TotalBases)), -10*log10(error/TotalBases)))
# qv_summary = data %>% group_by(ec,dataset) %>% filter(ec < 60) %>% summarise(numBases = sum(alnlen), qv = mean(qv))
g1 = ggplot(qv_summary) +
  # facet_grid(~dataset)+
  stat_smooth(aes(ec,qv,fill=dataset,col=dataset), alpha = 0.2) +
  geom_hline(yintercept = 30, lty = 3, alpha = 0.75) +
  geom_vline(xintercept = 15, lty = 3, alpha = 0.75) +
  geom_point(aes(ec,qv,size=numBases,col=dataset), alpha = 0.75) +
  coord_cartesian(ylim=c(15,50))+
  ylab("Mean Identity Phred")+
  xlab("Effective coverage")+
  theme_minimal()+
  scale_size_area(breaks=c(1 %o% 10^(1:20)))+
  theme(plot.title = element_text(hjust = 0.5, size=14), legend.position="right",
        axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"),
        legend.text = element_text(size=14), legend.title = element_text(size=14,face="bold"))+
  # theme(legend.position="bottom")+
  scale_x_continuous(limits=c(3,32), breaks=c(3,5,10,15,20,30))

ggsave("harmony.pdf",g1,width=25,height=15,dpi=200,units="cm",limitsize = FALSE)
