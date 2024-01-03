### this script was used to make scatter plots of piRNA genomic tile coverage and the table that summarises data from all libraries
library(ggplot2)
library(reshape2)

pdf.options(useDingbats = FALSE)

### setting the directory in R
setwd("${analysis}/tile-analysis")

### loading the table and dcast it
table=read.table("kb-tiles.txt", header=F)
colnames(table)=c("tile","value","cluster_or_library")
table_d <- dcast(table, tile~cluster_or_library, value.var="value")
table_d[is.na(table_d)] <- 0

### plotting tile coverage of esg-ts_GFP_ovaries_sRNA on X and esg-ts_GFP_aub-RNAi_ovaries_sRNA on Y.
### all other plots were made in the same way.
table_d2 <- subset(table_d, table_d$`esg-ts_GFP_ovaries_sRNA` > 1 & table_d$`esg-ts_GFP_aub-RNAi_ovaries_sRNA` > 1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`esg-ts_GFP_ovaries_sRNA`,
                 y=`esg-ts_GFP_aub-RNAi_ovaries_sRNA`),size=2,shape=1,alpha=0.8,colour="gray")+
  geom_point(data=subset(table_d2, Cluster42AB > 0.5),
             aes(x=`esg-ts_GFP_ovaries_sRNA`,
                 y=`esg-ts_GFP_aub-RNAi_ovaries_sRNA`),size=2,shape=1,alpha=0.6,colour="cyan")+
  geom_point(data=subset(table_d2, Cluster38CLeft > 0.5),
             aes(x=`esg-ts_GFP_ovaries_sRNA`,
                 y=`esg-ts_GFP_aub-RNAi_ovaries_sRNA`),size=2,shape=1,alpha=0.6,colour="orange")+
  geom_point(data=subset(table_d2, ClusterFlam > 0.5),
             aes(x=`esg-ts_GFP_ovaries_sRNA`,
                 y=`esg-ts_GFP_aub-RNAi_ovaries_sRNA`),size=2,shape=1,alpha=0.6,colour="green")+
  labs(title="ovary piRNAs vs gut ISCs piRNAs",
       x="esg-ts_GFP_ovaries_sRNA / RPKM", y="esg-ts_GFP_aub-RNAi_ovaries_sRNA / RPKM")+
  scale_x_log10(limits=c(1,20000))+
  scale_y_log10(limits=c(1,20000))+
  coord_fixed()+
  theme_bw()

pdf(file="kb-tiles.esg-ts_GFP_ovaries_sRNA.vs.esg-ts_GFP_aub-RNAi_ovaries_sRNA.pdf",width=8,height=8)
p
dev.off()

### printing the dcasted table
write.csv(table_d,file="kb-tiles.csv",row.names=F,col.names=T,quote=F)

### packages used in this script
# reshape2_1.4.4
# ggplot2_3.4.4
