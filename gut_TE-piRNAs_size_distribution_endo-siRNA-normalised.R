### size distribution and nucleotide composition of TE mapping reads normalised to endo-siRNAs abundance

#get in variables from bash
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

### use together with ${scripts}/gut_sRNAseq_analysis.sh
### reshape2_1.4.4 ggplot2_3.3.3 
library(reshape2)
library(ggplot2)

ENDO=as.numeric(ENDO)

### 1st minus --- START ---
v=c(DIRECTORY,"/",LIB,"/",LIB,"_TE_3MM_mappers.1st.minus.nucleotide.length.table")
vname=paste(v,collapse="")
table=read.table(vname,header=F)

colnames(table)=c("size","nucleotide","reads")

table$reads <- as.numeric(table$reads) / ENDO

p=c(DIRECTORY,"/",LIB,"/plots/",LIB,"_TE-mappers_1st_minus.endo-siRNA-normalised.pdf")
pname=paste(p,collapse="")

mains=c(LIB,", total TE mappers, 1st_minus")
mains=paste(mains,collapse="")

pdf(file=pname,width=8,height=6)
ggplot(table,aes(x=size,y=reads,fill=nucleotide))+
labs(title=mains, x="size", y="Reads / endo-siRNAs")+
geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

### printing out numbers in a csv table
t=c(DIRECTORY,"/",LIB,"/plots/",LIB,"_TE-mappers_1st_minus.endo-siRNA-normalised.csv")
tname=paste(t,collapse="")
table_d <- dcast(table,size~nucleotide,value.var="reads")
table_d[is.na(table_d)] <- 0
write.csv(table_d,file=tname,quote=F,row.names=F)

### 1st minus --- END ---


### 10th plus --- START ---
v=c(DIRECTORY,"/",LIB,"/",LIB,"_TE_3MM_mappers.10th.plus.nucleotide.length.table")
vname=paste(v,collapse="")
table=read.table(vname,header=F)

colnames(table)=c("size","nucleotide","reads")

table$reads <- as.numeric(table$reads) / ENDO

p=c(DIRECTORY,"/",LIB,"/plots/",LIB,"_TE-mappers_10th_plus.endo-siRNA-normalised.pdf")
pname=paste(p,collapse="")

mains=c(LIB,", total TE mappers, 10th_plus")
mains=paste(mains,collapse="")

pdf(file=pname,width=8,height=6)
ggplot(table,aes(x=size,y=reads,fill=nucleotide))+
labs(title=mains, x="size", y="Reads / endo-siRNAs")+
geom_bar(stat="identity")+
theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

### printing out numbers in a csv table
t=c(DIRECTORY,"/",LIB,"/plots/",LIB,"_TE-mappers_10th_plus.endo-siRNA-normalised.csv")
tname=paste(t,collapse="")
table_d <- dcast(table,size~nucleotide,value.var="reads")
table_d[is.na(table_d)] <- 0
write.csv(table_d,file=tname,quote=F,row.names=F)

### 10th plus --- END ---
