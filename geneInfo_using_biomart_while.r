library("biomaRt")
ensembl37 <- useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL') #select the ensembl version that contains hg19 (GRch37) 
ensembl <- useDataset("hsapiens_gene_ensembl",mart= ensembl37) #select human genes
attributes <- c("ensembl_gene_id", "ensembl_transcript_id","chromosome_name","percentage_gc_content","entrezgene","hgnc_symbol","hgnc_transcript_name","strand" ,"band","transcript_length") #select the attributes
filter <- "hgnc_symbol" #select the filter
gene_id <- read.table("gene_id.txt") #upload the values of our filter
gene_id <- as.matrix(gene_id)
geneInfo = data.frame() #creat an empty dataframe
#the followinf while loop to handel the gene_id arry 500 by 500 (500 based on the advise on ensmble website beacuse whene we call the getBM function by the whole gene_id array we will miss some values)
finished = FALSE
k=0
while(!finished)
{
	if((length(gene_id)-(k*500))>500)
	{
		gene_id_tmp = gene_id[c((500*k+1):(500*(k+1)))]
		k=k+1
	}
	else
	{
		gene_id_tmp = gene_id[c((500*k+1):length(gene_id))]
		finished = TRUE
	}
	geneInfo <- rbind(geneInfo,getBM(attributes = attributes, filters = filter, values = gene_id_tmp, mart = ensembl))
}
row.na <- apply(geneInfo, 1, function(x){is.na(x[4])}) #select the rows that contains NA in gc content column(also trancript_length)
geneInfo.f <- geneInfo[!row.na,] 
gene_list <- split(geneInfo.f, with(geneInfo.f, interaction(hgnc_symbol)), drop = TRUE) #split the data according the gene hgnc_symbolID
max_index<-lapply(gene_list,function(x){which.max(x[,10])}) #find the index of bigger transcript
#build a new data frame contains the required transcripts(the transcript with bigger length )
#I think we can enhance the following segmment of code
geneInfo.unique = data.frame()
for (i in 1:length(gene_list))
{
	new_row = gene_list[[i]][max_index[[i]],]
	geneInfo.unique <- rbind(geneInfo.unique, new_row) 
}
write.csv(geneInfo.unique, file = "geneInfo_bimomart.csv")
