### this script was used to make scatter plots of TE mRNA counts
library(ggplot2)
library(ggrepel)
library(reshape2)

pdf.options(useDingbats = FALSE)

### setting the directory in R
setwd("${analysis}")

### loading the table and dcast it
table=read.table("salmon_TE_counts.txt", header=F)
colnames(table)=c("TE","value","library")
table_d <- dcast(table, TE~library, value.var="value")

### label abundant TEs
high_in_ovary <- c("diver", "3S18", "blood", "Transpac", "rover", "Max-element", "Burdock", "I-element_new", "HMS-Beagle", "flea")
high_in_gut <- c("copia", "roo", "Dm88", "Doc")
high_in_both <- c("flea")

### plotting tile coverage of esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1 on X and esg-Gal4_ISCs-EBs_Pe_polyA_rep1 on Y.
### all other plots were made in the same way.
table_d2 <- subset(table_d, table_d$`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1` > 0.1 & table_d$`esg-Gal4_ISCs-EBs_Pe_polyA_rep1` > 0.1)
table_d2_o <- subset(table_d2, TE %in% high_in_ovary)
table_d2_g <- subset(table_d2, TE %in% high_in_gut)
table_d2_og <- subset(table_d2, TE %in% high_in_both)
rownames(table_d2_o) <- table_d2_o$TE
rownames(table_d2_g) <- table_d2_g$TE
rownames(table_d2_og) <- table_d2_og$TE

p<-ggplot(table_d2)+
  geom_point(aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                 y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),size=4,shape=1,alpha=0.8,colour="gray")+
  geom_point(data=table_d2_o,
             aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                 y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),size=4,shape=1,alpha=0.8,colour="blue")+
  geom_text_repel(data=table_d2_o,
                  aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                      y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),label=rownames(table_d2_o),colour="blue")+
  geom_point(data=table_d2_g,
             aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                 y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),size=4,shape=1,alpha=0.8,colour="red")+
  geom_text_repel(data=table_d2_g,
                  aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                      y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),label=rownames(table_d2_g),colour="red")+
  geom_point(data=table_d2_og,
             aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                 y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),size=4,shape=1,alpha=0.8,colour="green")+
  geom_text_repel(data=table_d2_og,
                  aes(x=`esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1`,
                      y=`esg-Gal4_ISCs-EBs_Pe_polyA_rep1`),label=rownames(table_d2_og),colour="green")+
  labs(title="esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1 vs esg-Gal4_ISCs-EBs_Pe_polyA_rep1, TE mRNAs",
       x="TPM / esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1", y="TPM / esg-Gal4_ISCs-EBs_Pe_polyA_rep1")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="TE-mRNA-scatter-plot.esg-Gal4_ISCs-EBs_Sucrose_polyA_rep1.vs.esg-Gal4_ISCs-EBs_Pe_polyA_rep1.pdf",width=6,height=6)
p
dev.off()

### printing the dcasted table
write.csv(table_d,file="TE-mRNA-counts.csv",row.names=F,col.names=T,quote=F)

### packages used in this script
# reshape2_1.4.4
# ggplot2_3.4.4
# ggrepel_0.9.4
