# load insulin analysis.
data <- read.table("merged_gene_rpkm_by_read.csv",header = T,sep="\t",stringsAsFactors = F)
geneinfo <- data[,1:5]
mat <- data[,-1:-5]

# raw counts
raw <- read.table("merged.geneCount.txt",header = T,sep="\t",stringsAsFactors = F)
raw <- raw[,-1:-5]
# filtering
# remove genes with group median raw counts in both group lower than 32
raw <- raw[,grep("Insulin|Control",colnames(raw),value = T ) ]

keep <- apply(raw,1,function(x) any(x>=32) )

rpkm.filtered <- mat[keep,grep("Insulin|Control",colnames(mat),value = T )]
geneinfo.filtered <- geneinfo[keep,]




rpkm.insulin.rawfiltered <- rpkm.filtered[, c("D3_Insulin_10min_rpkm","KK12_Insulin_20min_rpkm","D3_Insulin_40min_rpkm","D3_Insulin_80min_rpkm","KK11_Insulin_300min_rpkm") ]
rpkm.control.rawfiltered <- rpkm.filtered[,c("D3_Control_10min_rpkm","KK10_Control_20min_rpkm","D3_Control_40min_rpkm","D3_Control_80min_rpkm","KK13_Control_300min_rpkm")]
diff.rpkm.rawfilter <- log2(rpkm.insulin.rawfiltered+0.01) - log2(rpkm.control.rawfiltered+0.01)
write.table(data.frame(GeneID=geneinfo.filtered$GeneName,diff.rpkm.rawfilter),"/research/bsi/projects/staff_analysis/m092409/insulin/diff.expressed.filtered.raw.txt",quote=F,sep="\t",row.names = F)

~
