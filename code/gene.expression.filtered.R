# new stem analysis

# load insulin analysis. 
data <- read.table("insulin/merged_gene_rpkm_by_read.csv",header = T,sep="\t",stringsAsFactors = F)
geneinfo <- data[,1:7]
mat <- data[,-1:-7]

# raw counts 
raw <- read.table("insulin/merged_genecounts.txt",header = T,sep="\t",stringsAsFactors = F)
raw <- raw[,-1:-7]
# filtering 
# remove genes with group median raw counts < 32 in all samples 
raw <- raw[,grep("Insulin|Control",colnames(raw),value = T ) ]

keep <- apply(raw,1,function(x) any(x>=32) )  

rpkm.filtered <- mat[keep,grep("Insulin|Control",colnames(mat),value = T )]
geneinfo.filtered <- geneinfo[keep,]

medians.rpkm.filtered <- t(apply(rpkm.filtered,1,median.group,group))
keep2 <- apply(medians.rpkm.filtered,1, function(x) any(x>1))

rpkm.final.filtered <- data.matrix(rpkm.filtered[keep2,])
geneinfo.final.filtered <- geneinfo.filtered[keep2,]
rownames(rpkm.final.filtered) <- geneinfo.final.filtered$GeneName

rpkm.insulin <- rpkm.final.filtered[,grep("Insulin",colnames(rpkm.final.filtered),value = T)]
rpkm.insulin <- rpkm.final.filtered[, c("D3_Insulin_10min_rpkm","KK12_Insulin_20min_rpkm","D3_Insulin_40min_rpkm","D3_Insulin_80min_rpkm","KK11_Insulin_300min_rpkm") ]
rpkm.control <- rpkm.final.filtered[,c("D3_Control_10min_rpkm","KK10_Control_20min_rpkm","D3_Control_40min_rpkm","D3_Control_80min_rpkm","KK13_Control_300min_rpkm")]

diff.rpkm.change <- log2(rpkm.insulin+0.01) - log2(rpkm.control+0.01 )
head(diff.rpkm.change)
colnames(diff.rpkm.change)<-gsub("_Insulin|_rpkm","",colnames(diff.rpkm.change))
write.table(data.frame(GeneID=rownames(diff.rpkm.change),diff.rpkm.change),"/research/bsi/projects/staff_analysis/m092409/insulin/diff.expressed.filtered.rpkm.txt",quote=F,sep="\t",row.names = F)


rpkm.insulin.rawfiltered <- rpkm.filtered[, c("D3_Insulin_10min_rpkm","KK12_Insulin_20min_rpkm","D3_Insulin_40min_rpkm","D3_Insulin_80min_rpkm","KK11_Insulin_300min_rpkm") ]
rpkm.control.rawfiltered <- rpkm.filtered[,c("D3_Control_10min_rpkm","KK10_Control_20min_rpkm","D3_Control_40min_rpkm","D3_Control_80min_rpkm","KK13_Control_300min_rpkm")]
diff.rpkm.rawfilter <- log2(rpkm.insulin.rawfiltered+0.01) - log2(rpkm.control.rawfiltered+0.01)
write.table(data.frame(GeneID=geneinfo.filtered$GeneName,diff.rpkm.rawfilter),"/research/bsi/projects/staff_analysis/m092409/insulin/diff.expressed.filtered.raw.txt",quote=F,sep="\t",row.names = F)

cqn.diff <- read.table("/research/bsi/projects/staff_analysis/m092409/insulin/diff.allgene.cqn.rmT4T6.txt",header=T,sep="\t",stringsAsFactors = F,row.names = 1)
dim(cqn.diff)
cqn.diff.filter <-  cqn.diff[rownames(cqn.diff) %in% geneinfo.filtered$GeneName,]
write.table(data.frame(GeneID=rownames(cqn.diff.filter),cqn.diff.filter),"/research/bsi/projects/staff_analysis/m092409/insulin/diff.allgene.cqn.rmT4T6.filtered.txt",quote=F,sep="\t",row.names = F)



