### the code was used to analyse the small RNA sequencing data for gut and ovary samples from Drosophila melanogaster.

### software used
### fastx_toolkit
http://hannonlab.cshl.edu/fastx_toolkit/
### bowtie-1.2.3-linux-x86_64
https://github.com/BenLangmead/bowtie/releases/tag/v1.2.3
### samtools/1.10
https://github.com/samtools/samtools/releases/tag/1.10
### bedtools/2.28.0
https://github.com/arq5x/bedtools2/releases/tag/v2.28.0
### kent-ucsc tools
https://github.com/ucscGenomeBrowser/kent
### R/4.0.0
https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz
### weblogo 3.7.8
https://github.com/WebLogo/weblogo

### define the variables
scripts="<directory where the scripts are kept>"
references="<directory where the genome sequence file and other sequence files are kept>"
raw_data="<directory where the fastq files are kept>"
analysis="<directory where the output files are kept>"
lib_sRNA="<name of the small RNA seq library>"
### directory to save plots per library
mkdir -p ${analysis}/${lib_sRNA}/plots


### processing and analysing the small RNA sequencing libraries --- START ---

### below are the small RNA libraries processed in this script:
### from this study
w1118_midgut_Sucrose_sRNA_rep1
w1118_midgut_Sucrose_sRNA_rep2
w1118_midgut_Sucrose_sRNA_rep3
w1118_midgut_Pe_sRNA_rep1
w1118_midgut_Pe_sRNA_rep2
w1118_midgut_Pe_sRNA_rep3
aub-HN2-QC42_midgut_Sucrose_sRNA_rep1
aub-HN2-QC42_midgut_Sucrose_sRNA_rep2
aub-HN2-QC42_midgut_Sucrose_sRNA_rep3
aub-HN2-QC42_midgut_Pe_sRNA_rep1
aub-HN2-QC42_midgut_Pe_sRNA_rep2
aub-HN2-QC42_midgut_Pe_sRNA_rep3
esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep1
esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep2
esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep3
esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep1
esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep2
esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep3
esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep1
esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep2
esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep3
esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep1
esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep2
esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep3
esg-ts_GFP_ovaries_sRNA
esg-ts_GFP_aub-RNAi_ovaries_sRNA

### from Siudeja et al., 2021, EMBO J
### https://www.ebi.ac.uk/ena/browser/view/PRJEB41757
run_accession	sample_accession	experiment_title	fastq_ftp	fastq_aspera	sample_title
ERR4920944	SAMEA7690898	NextSeq 500 sequencing	ftp.sra.ebi.ac.uk/vol1/fastq/ERR492/004/ERR4920944/ERR4920944.fastq.gz	fasp.sra.ebi.ac.uk:/vol1/fastq/ERR492/004/ERR4920944/ERR4920944.fastq.gz	Gut-Pros>2xGFP-1
ERR4920945	SAMEA7690899	NextSeq 500 sequencing	ftp.sra.ebi.ac.uk/vol1/fastq/ERR492/005/ERR4920945/ERR4920945.fastq.gz	fasp.sra.ebi.ac.uk:/vol1/fastq/ERR492/005/ERR4920945/ERR4920945.fastq.gz	Gut-Pros>2xGFP-2

### preparing for the analyses, do this before processing individual libraries --- START ---
### step 0-1: generating a bowtie index for miscellaneous RNAs
### step 0-1-1: download miscRNA sequences from the FlyBase
mkdir -p ${references}/dm6/indices
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-miscRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-tRNA-r6.31.fasta.gz
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/releases/FB2019_06/dmel_r6.31/fasta/dmel-all-ncRNA-r6.31.fasta.gz

### step 0-1-2: making fasta including miRNA, rRNA, snRNA, snoRNA, and tRNA sequences
zcat ${references}/dm6/dmel-all-miscRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-miRNA-r6.31.fasta.gz ${references}/dm6/dmel-all-tRNA-r6.31.fasta.gz |\
fasta_formatter - -t | tr -d ";" |\
awk '{split($5,a,"="); if($2~"miRNA" || $2~"rRNA" || $2~"snRNA" || $2~"snoRNA" || $2~"tRNA") print ">"a[2]"\n"$NF
}' > ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta

### step 0-1-3: making a bowtie index
bowtie-build ${references}/dm6/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31.fasta ${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31


### step 0-2: generating a bowtie index for the genome
### step 0-2-1: download the genomic sequence of the dm6 r6.31 assembly of Drosophila melanogaster
wget --directory-prefix="${references}/dm6" http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.31_FB2019_06/fasta/dmel-all-chromosome-r6.31.fasta.gz
### only take the major chromosomes for the analysis
zcat ${references}/dm6/dmel-all-chromosome-r6.31.fasta.gz | fasta_formatter - -t |\
awk '{if($1=="2L" || $1=="2R" || $1=="3L" || $1=="3R" || $1=="4" || $1=="X" || $1=="Y" ) print ">chr"$1"\n"$NF}' > ${references}/dm6/dmel-all-chromosome-r6.31.nosmallChr.fasta

### step 0-2-2: make a bowtie index
bowtie-build ${references}/dm6/dmel-all-chromosome-r6.31.nosmallChr.fasta ${references}/dm6/indices/dmel-all-chromosome-r6.31.nosmallChr

### step 0-2-3: generate the size file
fasta_formatter -i ${references}/dm6/dmel-all-chromosome-r6.31.nosmallChr.fasta -t |\
awk '{print $1,length($NF)}' > ${references}/dm6/dmel-all-chromosome-r6.31.nosmallChr.sizes


### step 0-3: generating a bed file for excluding certain reads
### These include fragments of tRNAs, snoRNA and snRNA as well as siRNAs produced from the transgene that target aubergine.
### step 0-3-1: make a bed file of tRNA, snoRNA and snRNA annotations extended by 100 bp both upstream and downstream
zcat ${references}/dm6/dmel-all-tRNA-r6.31.fasta.gz | grep ">" | grep -v "mito" | awk 'BEGIN{ FS = ";"} {print $2}' | tr -d "loc=)" | tr ':(' ' ' |\
awk '{n=split($NF,a,"."); if($0 ~ "mpement") print "chr"$1,a[1]-101,a[n]+100,". . -";
else print "chr"$1,a[1]-101,a[n]+100,". . +"}' | sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/dm6/dmel-all-tRNA-r6.31.tRNA-extended.bed

### step 0-3-2: make a bed file of snoRNA and snRNA annotations extended by 100 bp both upstream and downstream
zcat ${references}/dm6/dmel-all-miscRNA-r6.31.fasta.gz | grep ">" | grep -e "snoRNA" -e "snRNA" | awk 'BEGIN{ FS = ";"} {print $2}' | tr -d "loc=)" | tr ':(' ' ' |\
awk '{n=split($NF,a,"."); if($0 ~ "mpement") print "chr"$1,a[1]-101,a[n]+100,". . -";
else print "chr"$1,a[1]-101,a[n]+100,". . +"}' | sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/dm6/dmel-all-tRNA-r6.31.snoRNA-snRNA-extended.bed

### step 0-3-3: make a bed file for the aubergine locus
printf "chr2L 10997819 11001476 . . -""\n""chr2L 10997819 11001476 . . +""\n" | tr ' ' '\t' > ${references}/dm6/aub.bed


### step 0-4: generating a bed tile for genome unique 0.5kb tiles from the dm6 genome assembly of Drosophila melanogaster
### step 0-4-1: obtain every 25mer from the genome, only sense
### convert lowercase to uppercase first
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${references}/dm6/dmel-all-chromosome-r6.31.nosmallChr.fasta |\
fasta_formatter - -t | awk '{for(i=0;i<=length($NF)-25;i++) print ">"$1"@"i"_+:"substr($NF,i+1,25)"\n"substr($NF,i+1,25)
}' > ${references}/dm6/Dmel_dm6_25mers.plus.fasta

### step 0-4-2: obtain genome unique mappers
index_genome="${references}/dm6/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 0 -m 1 -S ${index_genome} ${references}/dm6/Dmel_dm6_25mers.plus.fasta |\
samtools view -@ 16 -bS - | bamToBed -i - > ${references}/dm6/Dmel_dm6_25mers.plus.unique-mappers.bed

### step 0-4-3: count the unique mappers in 0.5kb bins
### take if the centre of 25mer sits between the coordinate of 1 to 500 as the first 0.5kb bin
### take as the second 0.5kb bin if it sits between 501 to 1000
awk '{split($4,a,"@"); split(a[2],b,"_"); TILE[$1"_"int((b[1]+12)/500)]++
} END {for(var in TILE) print var,TILE[var]}' ${references}/dm6/Dmel_dm6_25mers.plus.unique-mappers.bed > ${references}/dm6/Dmel_dm6_0.5kbtiles.counts

### step 0-4-4: make a bed file that contains all 0.5kb tiles that are gt 85% mappability
awk '{split($1,a,"_"); if($2>425) print a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_+ . +""\n"a[1],a[2]*500,(a[2]+1)*500,a[1]":"a[2]*500"-"(a[2]+1)*500"_- . -"
}' ${references}/dm6/Dmel_dm6_0.5kbtiles.counts |\
sort -k1,1 -k2,2n | tr ' ' '\t' > ${references}/dm6/Dmel_dm6_85percent_0.5kbtiles.25mers.bed
### ${references}/dm6/Dmel_dm6_85percent_0.5kbtiles.25mers.bed looks as follows:
### the fourth column breaks down to: <chromosome> @ <start> - <end> _ <strand> : <fraction of unique mappers>
chr2L   5500    6000    chr2L@5500-6000_+:0.964 1       +
chr2L   5500    6000    chr2L@5500-6000_+:0.964 1       -
chr2L   6000    6500    chr2L@6000-6500_+:1     1       +
...


### step 0-5: generating a bowtie index for endo siRNAs
zcat ${references}/dm6/dmel-all-ncRNA-r6.31.fasta.gz | fasta_formatter - -t | awk '{if($0 ~ "hpRNA") print ">"$1":hpRNA""\n"$NF}' > ${references}/dm6/dmel-hpRNA-r6.31.fasta
bowtie-build ${references}/dm6/dmel-hpRNA-r6.31.fasta ${references}/dm6/indices/dmel-hpRNA-r6.31


### step 0-6: make a bowtie index for the Drosophila melanogaster transposon sequences used in Senti G&D 2015, PMID: 26302790 (122 entries)
bowtie-build ${references}/TEs/122_dmel_TE_Senti2015.fa ${references}/TEs/indices/122_dmel_TE_Senti2015
fasta_formatter -i ${references}/TEs/122_dmel_TE_Senti2015.fa -t | awk '{print $1,length($2)}' > ${references}/TEs/122_dmel_TE_Senti2015.sizes


### common analyses for all small RNA libraries --- START ---

### step 1: trim the adapter sequence and collapse reads by sequence

### libraries from this study
### fastq_to_fasta, fasta_formatter, and fastx_clipper are from fastx_toolkit
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
# collapse reads of the same sequence to one entry while retaining the number of reads
fastx_clipper -c -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i - | fasta_formatter - -t |\
awk '{if(length($NF)>25 && length($NF)<49) READS[substr($NF,5,length($NF)-8)]++} END {
for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa

### libraries from Siudeja et al.
fastq_to_fasta -Q33 -i <(zcat ${raw_data}/${lib_sRNA}.fastq.gz) |\
# -c: discard non-clipped sequences, -l 18: discard sequences shorter than 18nt
fastx_clipper -c -a TGGAATTCTCGGGTGCCAAG -l 18 -i - | fasta_formatter - -t |\
awk '{READS[$NF]++} END {for(var in READS) print ">"var"@"READS[var]"\n"var}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa

### ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa looks as follows:
>TCCAGATGAACGGTTAAGTGTCCAAAAAG@12
TCCAGATGAACGGTTAAGTGTCCAAAAAG
>TACTTGAAAGAATCAGGGGCCAACCAT@5
TACTTGAAAGAATCAGGGGCCAACCAT
...


### step 2: map to the Drosophila melanogaster miscellaneous RNA (misc RNA) allowing up to 1MM
### run bowtie to map reads to the miscRNA and take unmapped reads
index_misc="${references}/dm6/indices/dmel-all-miRNA-rRNA-snRNA-snoRNA-tRNA-r6.31"
bowtie -f -v 1 -a -S ${index_misc} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed.fa \
--un ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-mapped.bed


### step 3: obtain genome unique mappers
### run bowtie to map reads to the genome, allowing up to 1MM, unique mappers
index_genome="${references}/dm6/indices/dmel-all-chromosome-r6.31.nosmallChr"
bowtie -f -v 1 -m 1 -S ${index_genome} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - |\
awk '{if($1=="chr2L" || $1=="chr2R" || $1=="chr3L" || $1=="chr3R" || $1=="chr4" || $1=="chrX") print}' |\
sort -k1,1 -k2,2n > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.bed


### step 4: exclude reads that intersect with tRNA, snoRNA, snRNA annotations as well as the aub locus
cat ${references}/dm6/dmel-all-tRNA-r6.31.tRNA-extended.bed ${references}/dm6/dmel-all-tRNA-r6.31.snoRNA-snRNA-extended.bed ${references}/dm6/aub.bed |\
bedtools intersect -v -s -f 0.5 -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.bed \
-b - > ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.tRNA-snoRNA-snRNA-aub-excluded.bed


### step 5: 0.5kb tile analysis for the genome unique piRNA (>22nt) mappers
### count the number of genome unique mappers that are longer than 22nt to be considered as piRNAs
READS=`(awk '{split($4,a,"@"); if(length(a[1])>22) count+=a[2]} END {print count
}' ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.tRNA-snoRNA-snRNA-aub-excluded.bed)`
### normalise counts per tile per one million genome unique mappers
bedtools intersect -wo -s -a ${analysis}/${lib_sRNA}/${lib_sRNA}_misc-unmapped_genome-unique-mappers.tRNA-snoRNA-snRNA-aub-excluded.bed \
-b ${references}/dm6/Dmel_dm6_85percent_0.5kbtiles.25mers.bed |\
awk -v READS=${READS} -v LIB=${lib_sRNA}  '{split($4,a,"@"); if(length(a[1])>22) TILE[$10]+=a[2]} END {for(var in TILE) print var,TILE[var]/READS*1000000,LIB
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}.tile-coverage.txt
### ${analysis}/${lib_sRNA}/${lib_sRNA}.tile-coverage.txt looks as follows:
### <tile coordinate> <tile coverage> <library name>


### step 6: map to the Drosophila melanogaster endogenous siRNAs allowing up to 1MM
index_endo_siRNAs="${references}/dm6/indices/dmel-hpRNA-r6.31"
bowtie -f -v 1 --all --best --strata -S ${index_endo_siRNAs} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -bS - | bamToBed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped_endo_siRNAs.bed


### step 7: mapping and characterisation of transposon small RNAs
### step 7-1: bowtie mapping to TEs, --all --best --strata option, allowing up to three mismatches
index_TE="${references}/TEs/indices/122_dmel_TE_Senti2015"
bowtie -f -v 3 --all --best --strata -S ${index_TE} ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped.fa |\
samtools view -b -S - | bedtools bamtobed -i - > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed

### step 7-2: counting the nucleotides occuring at the 1st and 10th positions of TE antisense (minus) and sense (plus) reads, respectively
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed |\
awk '{split($4,a,"@"); if($6=="-") FIRST[length(a[1])" "substr($4,1,1)]+=a[2]} END {for(var in FIRST) print var,FIRST[var]
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.1st.minus.nucleotide.length.table
awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed |\
awk '{split($4,a,"@"); if($6=="+") TENTH[length(a[1])" "substr($4,10,1)]+=a[2]} END {for(var in TENTH) print var,TENTH[var]
}' > ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.10th.plus.nucleotide.length.table

### step 7-3: normalising the counts by the number of reads mapping to endogenous siRNAs and making plots using R
ENDO=`(awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_collapsed_misc-unmapped_endo_siRNAs.bed | awk '{split($4,a,"@"); count+=a[2]} END {print count}')`
Rscript ${scripts}/gut_TE-piRNAs_size_distribution_endo-siRNA-normalised.R DIRECTORY=${analysis} LIB=${lib_sRNA} ENDO=${ENDO}

### step 7-4: measuring Uridine and Adenosine frequencies at the 1st and 10th positions of TE antisense and sense piRNAs using weblogo
mkdir -p ${analysis}/${lib_sRNA}/weblogo
TE_fasta="${references}/TEs/122_dmel_TE_Senti2015.fa"
cat ${references}/TEs/122_dmel_TE_Senti2015.sizes | while read TE length; do

  awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed | awk -v TE=${TE} '{if($1==TE && $3-$2>22) print}' |\
  awk '{split($4,a,"@"); split(a[2],b,":"); for(i=1;i<=b[1];i++) {if($6=="-") print $1,$3-6,$3+5,$4,$5,$6}}' |\
  awk -v LENGTH=${length} '{if($2>=0 && $3<=LENGTH) print}' | tr ' ' '\t' |\
  bedtools getfasta -s -fi ${TE_fasta} -tab -bed - | awk '{print ">"$1"\n"toupper($2)}' | tr 'T' 'U' >> ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_all-TE-AS-piRNAs.5end_11nt_window.fasta

  awk '!seen[$4]++' ${analysis}/${lib_sRNA}/${lib_sRNA}_TE_3MM_mappers.bed | awk -v TE=${TE} '{if($1==TE && $3-$2>22) print}' |\
  awk '{split($4,a,"@"); split(a[2],b,":"); for(i=1;i<=b[1];i++) {if($6=="+") print $1,$2+4,$2+15,$4,$5,$6}}' |\
  awk -v LENGTH=${length} '{if($2>=0 && $3<=LENGTH) print}' | tr ' ' '\t' |\
  bedtools getfasta -s -fi ${TE_fasta} -tab -bed - | awk '{print ">"$1"\n"toupper($2)}' | tr 'T' 'U' >> ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_all-TE-S-piRNAs.middle_11nt_window.fasta

done

for TYPE in all-TE-AS-piRNAs.5end_11nt_window all-TE-S-piRNAs.middle_11nt_window; do
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F pdf -n 50 -c classic -o ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.gt22_logo_prob.pdf
weblogo -U probability -A rna -f ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.fasta -F logodata -n 50 -c classic -s large -t ${lib}"_"${TE} -o ${analysis}/${lib_sRNA}/weblogo/${lib_sRNA}_${TYPE}.gt22_logo_prob.txt
done

### common analyses for all small RNA libraries --- END ---


### kb tile analysis, comparing between libraries --- START ---
### libraries that were used to generate plots
esg-ts_GFP_ovaries_sRNA # always used as a control wild type ovary library
esg-ts_GFP_aub-RNAi_ovaries_sRNA # aub-RNAi is expressed in the gut while libraries were made from ovarian RNA, effectively serving as a second wild type ovary library
### Below four libraries from sorted ISCs/EBs were compared against the wild type ovary library.
### One replicate each from two RNAi conditions, Sucrose and Pe feeding was used.
esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep2
esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep3
esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep1
esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep1

### coordiantes of known piRNA clusters in the dm6 gneome
### ${references}/dm6/Dmel_piRNAs-clusters.bed looks as follows:
chrX    21520428        21556793        Cluster20A      0       +
chrX    21631891        22139306        ClusterFlam     0       +
chr3L   23280548        23311746        Cluster80F      0       +
chr3L   23280548        23311746        Cluster80F      0       -
chr2R   6252727 6496564 Cluster42AB     0       +
chr2R   6252727 6496564 Cluster42AB     0       -
chr2L   20104254        20116386        Cluster38CLeft  0       +
chr2L   20104254        20116386        Cluster38CLeft  0       -
chr2L   20148381        20223595        Cluster38CRight 0       +
chr2L   20148381        20223595        Cluster38CRight 0       -

mkdir -p ${analysis}/tile-analysis
### step 8-1: collecting unique tiles in ClusterFlam (flamenco), Cluster42AB and Cluster38CLeft
cat ${references}/dm6/Dmel_piRNAs-clusters.bed | grep -e "ClusterFlam" -e "Cluster42AB" -e "Cluster38CLeft" |\
bedtools intersect -wo -s -a ${references}/dm6/Dmel_dm6_85percent_0.5kbtiles.25mers.bed -b - |\
awk '{print $4,$NF/500,$10}' > ${analysis}/tile-analysis/kb-tiles.txt 

### step 8-2: collecting tile coverage of the libraries of interests
for lib_sRNA in esg-ts_GFP_ovaries_sRNA esg-ts_GFP_aub-RNAi_ovaries_sRNA esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep2 esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep3 esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep1 esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep1; do
cat ${analysis}/${lib_sRNA}/${lib_sRNA}.tile-coverage.txt >> ${analysis}/tile-analysis/kb-tiles.txt
done

### step 8-3: making scatter plots of piRNA genomic tile coverage and the table that summarises data from all libraries
### ${analysis}/tile-analysis/kb-tiles.txt and ${scripts}/gut_piRNAs-genomic-tile-analysis.R was used.
### outputs: scatter plots used in the figures and the table ${analysis}/tile-analysis/kb-tiles.csv

### step 8-4: collecting normalised coverage for 8 most abundant tiles per piRNA cluster
CUTOFF_38C="130"
CUTOFF_42AB="320"
CUTOFF_Flam="400"
awk -v CUTOFF_38C=${CUTOFF_38C} -v CUTOFF_42AB=${CUTOFF_42AB} -v CUTOFF_Flam=${CUTOFF_Flam} 'BEGIN{FS=","; OFS=","; print "tile_coordinate,Cluster38CLeft,Cluster42AB,ClusterFlam,esg-ts_GFP_aub-RNAi_ovaries_sRNA,esg-ts_GFP_ISCs-EBs_Sucrose_sRNA_rep2,esg-ts_GFP_ISCs-EBs_Pe_sRNA_rep3,esg-ts_GFP_aub-RNAi_ISCs-EBs_Sucrose_sRNA_rep1,esg-ts_GFP_aub-RNAi_ISCs-EBs_Pe_sRNA_rep1"
} {if(NR>1 && $5>CUTOFF_38C && $6>CUTOFF_38C && $7>CUTOFF_38C && $8>CUTOFF_38C && $9>CUTOFF_38C && $10>CUTOFF_38C && $2>0.5) print $1,$2,$3,$4,$6/$5,$7/$5,$8/$5,$9/$5,$10/$5;
else if(NR>1 && $5>CUTOFF_42AB && $6>CUTOFF_42AB && $7>CUTOFF_42AB && $8>CUTOFF_42AB && $9>CUTOFF_42AB && $10>CUTOFF_42AB && $3>0.5) print $1,$2,$3,$4,$6/$5,$7/$5,$8/$5,$9/$5,$10/$5;
else if(NR>1 && $5>CUTOFF_Flam && $6>CUTOFF_Flam && $7>CUTOFF_Flam && $8>CUTOFF_Flam && $9>CUTOFF_Flam && $10>CUTOFF_Flam && $4>0.5) print $1,$2,$3,$4,$6/$5,$7/$5,$8/$5,$9/$5,$10/$5
}' ${analysis}/tile-analysis/kb-tiles.csv > ${analysis}/tile-analysis/kb-tiles.abundant-cluster-tiles.csv

### processing and analysing the small RNA sequencing libraries --- END ---
