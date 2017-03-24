library(tm)
library(ggplot2)
library(gridExtra)
library(plyr)
library(shiny)
library(reshape2)
library(reshape)
library(stringr)

setwd("/Users/mgehan/Documents/setaria/setaria-circadian/diel-explorer-112016/")

ldhhf<-read.table(file='data/JTK.LDHHF.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
llhcf<-read.table(file='data/JTK.LLHCF.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)
annotation<-read.table(file='data/setaria-annotation.txt',sep='\t',header=TRUE, stringsAsFactors = FALSE,strip.white = TRUE)

ldhhf$dataset<-"ldhhf"
llhcf$dataset<-"llhcf"

ldhhf.llhcf<-rbind(ldhhf,llhcf)

ldhhf.llhcf$species<-"setaria.viridis"

circadian.annotation<-merge(ldhhf.llhcf,annotation,by.x="GENEID",by.y="transcriptName", all.x=TRUE)
ldhhf.llhcf<-circadian.annotation
ldhhf.llhcf$ortholog<-NA
ldhhf.llhcf$ortholog<-paste(ldhhf.llhcf$Best.hit.arabi.name,ldhhf.llhcf$Best.hit.rice.name,sep=",")

ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==0]<-"00"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==1]<-"01"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==2]<-"02"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==3]<-"03"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==4]<-"04"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==5]<-"05"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==6]<-"06"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==7]<-"07"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==8]<-"08"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==9]<-"09"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==10]<-"10"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==11]<-"11"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==12]<-"12"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==13]<-"13"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==14]<-"14"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==15]<-"15"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==16]<-"16"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==17]<-"17"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==18]<-"18"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==19]<-"19"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==20]<-"20"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==21]<-"21"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==22]<-"22"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==23]<-"23"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==24]<-"24"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==25]<-"25"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==26]<-"26"
ldhhf.llhcf$LAG[ldhhf.llhcf$LAG==27]<-"27"

ldhhf.llhcf$bhq.cutoff<-NA
ldhhf.llhcf$bhq.cutoff[ldhhf.llhcf$BH.Q<=1]<- "1e+00"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.1)]<-"1e-01"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.01)]<-"1e-02"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.001)]<-"1e-03"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.0001)]<-"1e-04"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.00001)]<-"1e-05"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.000001)]<-"1e-06"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.0000001)]<-"1e-07"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.00000001)]<-"1e-08"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.000000001)]<-"1e-09"
ldhhf.llhcf$bhq.cutoff[(ldhhf.llhcf$BH.Q<=0.0000000001)]<-"1e-10"

ldhhf.llhcf$adjp.cutoff<-NA
ldhhf.llhcf$adjp.cutoff[ldhhf.llhcf$ADJ.P<=1]<-"1e+00"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.1)]<-"1e-01"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.01)]<-"1e-02"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.001)]<-"1e-03"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.0001)]<-"1e-4"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.00001)]<-"1e-5"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.000001)]<-"1e-6"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.0000001)]<-"1e-7"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.00000001)]<-"1e-8"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.000000001)]<-"1e-9"
ldhhf.llhcf$adjp.cutoff[(ldhhf.llhcf$ADJ.P<=0.0000000001)]<-"1e-10"

ldhhf.llhcf.max1<-subset(ldhhf.llhcf,select=c(ZT02_rep1,ZT04_rep1,ZT06_rep1,ZT08_rep1,ZT10_rep1,ZT12_rep1,
                                              ZT14_rep1,ZT16_rep1,ZT18_rep1,ZT20_rep1,ZT22_rep1,ZT24_rep1,ZT26_rep1,
                                              ZT28_rep1,ZT30_rep1,ZT32_rep1,ZT34_rep1,ZT36_rep1,ZT38_rep1,ZT40_rep1,
                                              ZT42_rep1,ZT44_rep1,ZT46_rep1,ZT48_rep1))
ldhhf.llhcf.max2<-subset(ldhhf.llhcf,select=c(ZT02_rep2,ZT04_rep2,ZT06_rep2,ZT08_rep2,ZT10_rep2,ZT12_rep2,
                                              ZT14_rep2,ZT16_rep2,ZT18_rep2,ZT20_rep2,ZT22_rep2,ZT24_rep2,ZT26_rep2,
                                              ZT28_rep2,ZT30_rep2,ZT32_rep2,ZT34_rep2,ZT36_rep2,ZT38_rep2,ZT40_rep2,
                                              ZT42_rep2,ZT44_rep2,ZT46_rep2,ZT48_rep2))

ldhhf.llhcf$max.rep1<-NA
ldhhf.llhcf$max.rep1<-apply(ldhhf.llhcf.max1,1,max)

ldhhf.llhcf$max.rep2<-NA
ldhhf.llhcf$max.rep2<-apply(ldhhf.llhcf.max2,1,max)

write.table(ldhhf.llhcf, file="data/setaria-ldhhf.llhcf.txt",quote=FALSE,sep='\t')

ldhhf.llhcf.norm<-subset(ldhhf.llhcf, select=-c(ZT02_rep1,ZT04_rep1,ZT06_rep1,ZT08_rep1,ZT10_rep1,ZT12_rep1,
                                        ZT14_rep1,ZT16_rep1,ZT18_rep1,ZT20_rep1,ZT22_rep1,ZT24_rep1,ZT26_rep1,
                                        ZT28_rep1,ZT30_rep1,ZT32_rep1,ZT34_rep1,ZT36_rep1,ZT38_rep1,ZT40_rep1,
                                        ZT42_rep1,ZT44_rep1,ZT46_rep1,ZT48_rep1,ZT02_rep2,ZT04_rep2,ZT06_rep2,ZT08_rep2,ZT10_rep2,ZT12_rep2,
                                        ZT14_rep2,ZT16_rep2,ZT18_rep2,ZT20_rep2,ZT22_rep2,ZT24_rep2,ZT26_rep2,
                                        ZT28_rep2,ZT30_rep2,ZT32_rep2,ZT34_rep2,ZT36_rep2,ZT38_rep2,ZT40_rep2,
                                        ZT42_rep2,ZT44_rep2,ZT46_rep2,ZT48_rep2))

ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT02_rep1= (ldhhf.llhcf$ZT02_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT04_rep1= (ldhhf.llhcf$ZT04_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT06_rep1= (ldhhf.llhcf$ZT06_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT08_rep1= (ldhhf.llhcf$ZT08_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT10_rep1= (ldhhf.llhcf$ZT10_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT12_rep1= (ldhhf.llhcf$ZT12_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT14_rep1= (ldhhf.llhcf$ZT14_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT16_rep1= (ldhhf.llhcf$ZT16_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT18_rep1= (ldhhf.llhcf$ZT18_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT20_rep1= (ldhhf.llhcf$ZT20_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT22_rep1= (ldhhf.llhcf$ZT22_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT24_rep1= (ldhhf.llhcf$ZT24_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT26_rep1= (ldhhf.llhcf$ZT26_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT28_rep1= (ldhhf.llhcf$ZT28_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT30_rep1= (ldhhf.llhcf$ZT30_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT32_rep1= (ldhhf.llhcf$ZT32_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT34_rep1= (ldhhf.llhcf$ZT34_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT36_rep1= (ldhhf.llhcf$ZT36_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT38_rep1= (ldhhf.llhcf$ZT38_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT40_rep1= (ldhhf.llhcf$ZT40_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT42_rep1= (ldhhf.llhcf$ZT42_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT44_rep1= (ldhhf.llhcf$ZT44_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT46_rep1= (ldhhf.llhcf$ZT46_rep1)/(ldhhf.llhcf$max.rep1))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT48_rep1= (ldhhf.llhcf$ZT48_rep1)/(ldhhf.llhcf$max.rep1))

ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT02_rep2= (ldhhf.llhcf$ZT02_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT04_rep2= (ldhhf.llhcf$ZT04_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT06_rep2= (ldhhf.llhcf$ZT06_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT08_rep2= (ldhhf.llhcf$ZT08_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT10_rep2= (ldhhf.llhcf$ZT10_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT12_rep2= (ldhhf.llhcf$ZT12_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT14_rep2= (ldhhf.llhcf$ZT14_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT16_rep2= (ldhhf.llhcf$ZT16_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT18_rep2= (ldhhf.llhcf$ZT18_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT20_rep2= (ldhhf.llhcf$ZT20_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT22_rep2= (ldhhf.llhcf$ZT22_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT24_rep2= (ldhhf.llhcf$ZT24_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT26_rep2= (ldhhf.llhcf$ZT26_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT28_rep2= (ldhhf.llhcf$ZT28_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT30_rep2= (ldhhf.llhcf$ZT30_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT32_rep2= (ldhhf.llhcf$ZT32_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT34_rep2= (ldhhf.llhcf$ZT34_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT36_rep2= (ldhhf.llhcf$ZT36_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT38_rep2= (ldhhf.llhcf$ZT38_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT40_rep2= (ldhhf.llhcf$ZT40_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT42_rep2= (ldhhf.llhcf$ZT42_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT44_rep2= (ldhhf.llhcf$ZT44_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT46_rep2= (ldhhf.llhcf$ZT46_rep2)/(ldhhf.llhcf$max.rep2))
ldhhf.llhcf.norm<- transform(ldhhf.llhcf.norm, ZT48_rep2= (ldhhf.llhcf$ZT48_rep2)/(ldhhf.llhcf$max.rep2))

write.table(ldhhf.llhcf.norm, file="data/setaria-ldhhf.llhcf.norm.txt",quote=FALSE,sep='\t')

