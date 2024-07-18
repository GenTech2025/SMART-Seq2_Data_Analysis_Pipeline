#In this test code we use scSEGIndex to fit model to a dataset and then compare 
#the top 1000 stable genes between two single celldata sets- one from ES (from class) one 
#from early development.  We compare the rank change between the top1000 SEG genes between
#the two samples and rank genes by the amount and direction of change between studies
#This rank is then used by fgsea as ranked input and compared to the molecular signatures db.
#The list also carries out rank tests to permute the sorted list and search for enrichment
#We then make some test data for David

#Finally the scSEGIndex is used to fit two different types of ES cells and the 
#enrichment comparison repeated!

#This code is testcode to show proof-of-concept.  It needs some extra organisation and 
#tidying before being considered "production" code.

#Data is from ES study used in the FGT class- (E-MTAB-2600)

#counts table from E4.5, E5.5, E6.5, E7.5) sequenced with scNMT-seq from Argelaguet et al. Nature. 2019 (GSE121708)

#scRNA data was not re-processed from Fastq for this test-used FGT class alignment & the counts from the GSE121708 record 

#data files used are on bioinfmsc6, in the /shared_files FGT_T9 folder for E-MTAB-2600, GSE121708 for GSE121708

#Simon Tomlinson May 2024



library(SingleCellExperiment)
library(scMerge)
#Data from E4.5, E5.5, E6.5, E7.5) sequenced with scNMT-seq from Argelaguet et al. Nature. 2019
counts <-read.table("../TestData/counts.txt",sep="\t",header=T,row.names=1)
tab_annotation <-read.table("../TestData/sample_metadata.txt", header=1, sep = "\t" )
counts_ordered <-counts[,tab_annotation$id_rna]

colSums(counts_ordered)

scedev <-SingleCellExperiment(counts_ordered)

#simple normalisation based on size factors
libsizes <- colSums(counts_ordered)
size.factors <- libsizes/mean(libsizes)
logcounts(scedev) <- log2(t(t(counts_ordered)/size.factors) + 1)
assayNames(scedev)

#check for NA values
sum(is.na(rowSums(logcounts(scedev))))
sum(is.na(colSums(logcounts(scedev))))
which(is.na(colSums(logcounts(scedev))))

annot_factors <-tab_annotation$stage
names(annot_factors) <-tab_annotation$id_rna

colsums <-colSums(counts_ordered)

tiff("../Plots/Outline_Code/colsums.tiff", width = 6, height = 6, units = "in", res = 500)
hist(colsums, breaks = 40)
dev.off()

#filter log counts from samples with over 100k reads
logcounts_filtered <- logcounts(scedev)[,colsums>100000]

tiff("../Plots/Outline_Code/log_filtered.tiff", width = 6, height = 6, units = "in", res = 500)
hist(colSums(logcounts_filtered))
dev.off()

sum(is.na(logcounts_filtered))
#OK to proceed

#adjust annotation so it matches the filtered data
annot_factors_filtered <-annot_factors[(names(annot_factors) %in% colnames(logcounts_filtered))]

#still need to cut out small libraries
res <-scSEGIndex(logcounts_filtered, cell_type = annot_factors_filtered)
res_ordered <-res[order(res$segIdx,decreasing = TRUE),]

#top match
tiff("../Plots/Outline_Code/topmatch.tiff", width = 6, height = 6, units = "in", res = 500)
hist(logcounts_filtered["ENSMUSG00000069744",],breaks = 40)
dev.off()

#bottom match
tiff("../Plots/Outline_Code/bottommatch.tiff", width = 6, height = 6, units = "in", res = 500)
hist(logcounts_filtered["ENSMUSG00000056300",],breaks = 40)
dev.off()

#do process with test data set to see how the stable genes differ between studies (this was later done in R below)

sce_escounts <- exprs(logNormCounts(sce))
sum(is.na(sce_escounts))

slotNames(sce)

#OK run scSEGIndex to find invariant genes
res_ES <-scSEGIndex(sce_escounts, cell_type = sce$Type)
res_ES_ordered <-res_ES[order(res_ES$segIdx,decreasing = TRUE),]

#top match
tiff("../Plots/Outline_Code/old_top.tiff", width = 6, height = 6, units = "in", res = 500)
hist(sce_escounts["ENSMUSG00000024991",],breaks = 40)
dev.off()

#bottom match
tiff("../Plots/Outline_Code/old_bottommatch.tiff", width = 6, height = 6, units = "in", res = 500)
hist(sce_escounts["ENSMUSG00000025839",],breaks = 40)
dev.off()

#write out the two tables for further analysis
write.csv(res_ES_ordered,file="../OutputFiles/res_ES_ordered.csv")
write.csv(res_ordered,file = "../OutputFiles/res_ordered_dev.csv")


#some preprocessing and then load back the results
es_rank <-read.csv("../Data/Outline_Data/extra-files/es_rank.csv",row.names =1)
dev_rank <-read.csv("../Data/Outline_Data/extra-files/dev_rank.csv", row.names =1)

es <-(es_rank)[,1]
names(es) <-rownames(es_rank)

dev <-(dev_rank)[,1]
names(dev) <-rownames(dev_rank)
dev <-dev[dev!="#N/A"]
dev2 <-as.numeric(dev)
names(dev2) <-names(dev)
dev <- dev2

png(filename = "../Plots/Outline_Code/es_barplot.png")
barplot(sort(es, decreasing = TRUE), main = "Barplot of Sorted ES", col = "blue")
dev.off()

png(filename = "../Plots/Outline_Code/dev_barplot.png")
barplot(sort(dev, decreasing = TRUE), main = "Barplot of Sorted DEV", col = "red")
dev.off()

####################################################################### RAN TILL HERE #########################

#perform enrichment
library(fgsea)

#example pathways loaded from R libraries (to test how to run the script)
data(examplePathways)
data(exampleRanks)
fgseaRes <- fgsea(examplePathways, exampleRanks, maxSize=500) #this runs the enrichment analysis

#BiocManager::install("msigdb")
library("msigdb")
library(ExperimentHub)
library(GSEABase)

eh = ExperimentHub()
query(eh , 'msigdb')

msigdb.mm = getMsigdb(org = 'mm', id = 'SYM', version = '7.5.1') #grab latest mouse data
msigdb.mm = appendKEGG(msigdb.mm) #add KEGG
msigdb.mm.list <- geneIds(msigdb.mm)

#fgsea using full signatures...
fgsea_es <- fgsea(msigdb.mm.list, es, maxSize=500,nPermSimple = 100000)

fgsea_dev <- fgsea(msigdb.mm.list, dev, maxSize=500,nPermSimple = 100000)

#fgsea using full msigdb... maxSize =500 may block chip enrichment matches (max list size often >500 in that case)
fgsea_es <- fgsea(msigdb.mm.list, es, nPermSimple = 100000)

fgsea_dev <- fgsea(msigdb.mm.list, dev, nPermSimple = 100000)

fgsea_es_ordered <-fgsea_es[order(padj, -abs(NES)), ]
fgsea_dev_ordered <-fgsea_dev[order(padj, -abs(NES)), ]

#Review of results showed some interesting signals
#https://www.cell.com/cell-metabolism/pdf/S1550-4131(13)00249-0.pdf -track across cell types
#https://doi.org/10.1016/j.cmet.2018.11.007 -mito dynamics...

fgsea_es_ordered[1:20,]$leadingEdge
fgsea_dev_ordered[1:20,]$leadingEdge

#eg fgsea_dev Rpl23, Rpl14, Rpl23a- paper Ribosomal proteins regulate 2-cell-stage transcriptome in mouse embryonic stem cells
#eg fgsea_dev Rpl23, Rpl14, Rpl23a Ribosomal proteins regulate 2-cell-stage transcriptome in mouse embryonic stem cells



listSubCollections(msigdb.mm) #see what other data we can find

#but we need to permute the gene set ranks to ensure we are not seeing biased results from picking
#the top ranked stable genes- any permute of these might give the same result?

#make randomised es gene set
es_rand <-es
names(es_rand) <-NULL
es_rand <-sample(es_rand)
names(es_rand) <-names(es)

#make randomised dev
dev_rand <-dev
names(dev_rand) <-NULL
dev_rand <-sample(dev_rand)
names(dev_rand) <-names(dev)


#run enrichment on these randomised sets -best to loop and summarise many results- here just a test
fgsea_es_rand <- fgsea(msigdb.mm.list, es_rand, maxSize=500,nPermSimple = 100000)
head(fgsea_es_rand[order(padj, -abs(NES)), ], n=20)


#run enrichment on these randomised sets
fgsea_dev_rand <- fgsea(msigdb.mm.list, dev_rand, maxSize=500,nPermSimple = 100000)
head(fgsea_dev_rand[order(padj, -abs(NES)), ], n=20)


#test of DAVID random 1000 element lists and fgsea 
test_names <- sample(rownames(sce_escounts))
test_names1000 <- test_names[1:1000]
test_names[grep("Notch", test_names, ignore.case = TRUE)]
test_names1000_Notch1 <-c(test_names1000,"Notch1")

#These lists go into DAVID and run below in fgsea

#prepare lists fgsea
rnames <-1:length(test_names1000)
rvalues <- test_names1000
test_names1000r <-rnames
names(test_names1000r)<-rvalues

rnames <-1:length(test_names1000_Notch1)
rvalues <- test_names1000_Notch1
test_names1000_Notch1r <- rnames
names(test_names1000_Notch1r) <- rvalues

#confirm notch is there
test_names1000_Notch1[grep("Notch", test_names1000_Notch1, ignore.case = TRUE)]

fgsea_test_names1000r <- fgsea(msigdb.mm.list, test_names1000r, maxSize=500,nPermSimple = 100000)
head(fgsea_test_names1000r[order(padj, -abs(NES)), ], n=20)

fgsea_test_names1000_Notch1r <- fgsea(msigdb.mm.list, test_names1000_Notch1r, maxSize=500,nPermSimple = 100000)
head(fgsea_test_names1000_Notch1r[order(padj, -abs(NES)), ], n=20)
#leading edge is Notch1

#check leading edge analysis- leading edge are core genes from pathways and signatures
fgsea_test_names1000_Notch1r_ordered <-fgsea_test_names1000_Notch1r[order(padj, -abs(NES)), ]
fgsea_test_names1000_Notch1r_ordered[1:5,]$leadingEdge


write.csv(test_names1000,file="test_names1000.csv")
write.csv(test_names, file = "test_names.csv")
write.csv(test_names1000_Notch1,file="test_names1000_Notch.csv")

plotEnrichment(pathwaysH[["<some Pathway"]], ranks) # needs to be plotted

#try fitting model within the scES data from class and cf serum and 2i (Types of cells)
sce$Type

sce_serum <-sce[,sce$Type =="serum"]
sce_2i <-sce[,sce$Type =="2i"]

#grab the log counts
sce_escounts_serum <- exprs(logNormCounts(sce_serum))
sce_escounts_2i <- exprs(logNormCounts(sce_2i))

#Input must be sorted this way or otherwise  gene names get scrambled!!
sce_escounts_serum <- sce_escounts_serum[order(rownames(sce_escounts_serum),decreasing =FALSE),]
sce_escounts_2i <- sce_escounts_2i[order(rownames(sce_escounts_2i),decreasing =FALSE),]

#test find genes changing stability
res_serum <-scSEGIndex(sce_escounts_serum, return_all = TRUE)
res_serum_ordered <-res_serum[order(res_serum$segIdx,decreasing = FALSE),]

res_2i <-scSEGIndex(sce_escounts_2i, return_all = TRUE)
res_2i_ordered <-res_2i[order(res_2i$segIdx,decreasing = FALSE),]

#add rank to each table
#add rank in other cell type
#store the rank difference
#sort by rank difference

#Sumo ppotency https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8604812/ - in some lists not others...

#add a rank to the table
res_serum_ordered$segIdx_serum <-rank(-1 * res_serum_ordered$segIdx,ties.method = "average")
res_2i_ordered$segIdx_2i <- rank(-1 * res_2i_ordered$segIdx, ties.method = "average")

#add cross element- put 2i and serum results in both tables
res_serum_ordered$segIdx_2i <- res_2i_ordered[rownames(res_serum_ordered),]$segIdx_2i
res_2i_ordered$segIdx_serum <- res_serum_ordered[rownames(res_2i_ordered),]$segIdx_serum

#calculate difference- 2i serum ranks each table
res_serum_ordered$segIdx_serumv2i <- res_serum_ordered$segIdx_serum - res_serum_ordered$segIdx_2i
res_2i_ordered$segIdx_2ivserum <- res_2i_ordered$segIdx_2i - res_2i_ordered$segIdx_serum

#tidy up mistakes!
#res_serum_ordered <-res_serum_ordered[-12]

#order again by segIdx- the stability rank
res_2i_ordered <-res_2i_ordered[order(res_2i_ordered$segIdx,decreasing = TRUE),]
res_serum_ordered <-res_serum_ordered[order(res_serum_ordered$segIdx,decreasing = TRUE),]

#take the top 1000, then order difference between 2i and serum
res_2i_ordered1000 <-res_2i_ordered[1:1000,]
res_serum_ordered1000 <-res_serum_ordered[1:1000,]

res_2i_ordered1000 <- res_2i_ordered1000[order(res_2i_ordered1000$segIdx_2ivserum, decreasing = TRUE),]
res_serum_ordered1000 <- res_serum_ordered1000[order(res_serum_ordered1000$segIdx_serumv2i, decreasing = TRUE),]

#make the input objects for enrichment- needs vector with rank values, names genes as row labels
res_2i_ordered1000_lst <-res_2i_ordered1000$segIdx_2ivserum
names(res_2i_ordered1000_lst) <- rownames(res_2i_ordered1000)

res_serum_ordered1000_lst <-res_serum_ordered1000$segIdx_serumv2i
names(res_serum_ordered1000_lst) <- rownames(res_serum_ordered1000)

#write the list
write.csv(res_2i_ordered1000_lst,file="res_2i_ordered1000_lst.csv")
write.csv(res_serum_ordered1000_lst,file="res_serum_ordered1000_lst.csv")

#run the enrichment analysis
fgsea_2i_ordered1000_lstr <- fgsea(msigdb.mm.list, res_2i_ordered1000_lst, nPermSimple = 100000)
fgsea_2i_ordered1000_lstr <- fgsea_2i_ordered1000_lstr[order(padj, -abs(NES)), ]
fgsea_2i_ordered1000_lstr[1:20]
fgsea_2i_ordered1000_lstr[1:20,]$leadingEdge

fgsea_serum_ordered1000_lstr <- fgsea(msigdb.mm.list, res_serum_ordered1000_lst, nPermSimple = 100000)
fgsea_serum_ordered1000_lstr <- fgsea_serum_ordered1000_lstr[order(padj, -abs(NES)), ]
fgsea_serum_ordered1000_lstr[1:20,]
fgsea_serum_ordered1000_lstr[1:20,]$leadingEdge

#finished....

sessionInfo()

