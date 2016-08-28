gInfo= read.csv("geneInfo_bimomart.csv",header=T,row.names=1,stringsAsFactors=FALSE)
colnames(gInfo)= c("ensembl_gene_id" ,"ensembl_transcript_id", "chr","gcContent","entrezgene","hgnc_symbol","hgnc_transcript_name","strand","band","geneLength")
geneInfo = gInfo[,c(10,4,3,1,2,5,6,7,8,9)]
save(gInfo, file = "geneInfo_biomart.rda")
