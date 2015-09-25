library('edgeR')

gly7a <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_gly7a.counts", row.names=1)
gly7b <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_gly7b.counts", row.names=1)
gly5a <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_gly5a.counts", row.names=1)
gly5b <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_gly5b.counts", row.names=1)
pyr7a <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_pyr7a.counts", row.names=1)
pyr7b <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_pyr7b.counts", row.names=1)
pyr5a <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_pyr5a.counts", row.names=1)
pyr5b <- read.table("~/TranscriptomicsWorkshop/TranscriptAbund/abundcount_pyr5b.counts", row.names=1)

colnames(gly7a) <- "Glycerol7a"
colnames(gly7b) <- "Glycerol7b"
colnames(gly5a) <- "Glycerol5a"
colnames(gly5b) <- "Glycerol5b"
colnames(pyr7a) <- "Pyruvate7a"
colnames(pyr7b) <- "Pyruvate7b"
colnames(pyr5a) <- "Pyruvate5a"
colnames(pyr5b) <- "Pyruvate5b"

gly7_vs_gly5 <- cbind(gly7a, gly7b, gly5a, gly5b) #compare glycerol pH 7.0 and pH 5.7
gly7_vs_pyr5 <- cbind(gly7a, gly7b, pyr5a, pyr5b) #compare glycerol pH 7.0 and pyruvate pH 5.7
pyr7_vs_pyr5 <- cbind(pyr7a, pyr7b, pyr5a, pyr5b) #compare pyruvate pH 7.0 and 5.7
gly7_vs_pyr7 <- cbind(gly7a, gly7b, pyr7a, pyr7b) #comparing the two conditions with which TB grows normally (shouldn't really be many differences)

group_ref_vs_trt <- c(1,1,2,2)

g7_vs_g5 <- DGEList(counts=gly7_vs_gly5, group=group_ref_vs_trt) #create DGEList data frame comparing glycerol pH 7 vs pH 5.7
g7_vs_p5 <- DGEList(counts=gly7_vs_pyr5, group=group_ref_vs_trt) #create DGEList data frame comparing glycerol pH 7 vs pyruvate pH 5.7
p7_vs_p5 <- DGEList(counts=pyr7_vs_pyr5, group=group_ref_vs_trt) #create DGEList data frame comparing pyruvate pH 7 vs pH 5.7
g7_vs_p7 <- DGEList(counts=gly7_vs_pyr7, group=group_ref_vs_trt) #create DGEList data frame comparing glycerol pH 7 vs pyruvate pH 7

prefilter_g7_vs_g5 <- rowSums(cpm(g7_vs_g5)>5) >= 2 #create the filter for glycerol pH 7 vs pH 5.7
prefilter_g7_vs_p5 <- rowSums(cpm(g7_vs_p5)>5) >= 2 #create the filter for glycerol pH 7 vs pyruvate pH 5.7
prefilter_p7_vs_p5 <- rowSums(cpm(p7_vs_p5)>5) >= 2 #create the filter for pyruvate pH 7 vs pyruvate pH 5.7
prefilter_g7_vs_p7 <- rowSums(cpm(g7_vs_p7)>5) >= 2 #create the filter for glycerol pH 7 vs pyruvate pH 7

#Apply the filter to the data
g7_vs_g5 <- g7_vs_g5[prefilter_g7_vs_g5,] #applying our filter we created earlier (its a vector that corresponds to which genes meet our filter criteria)
g7_vs_p5 <- g7_vs_p5[prefilter_g7_vs_p5,] #applying our filter we created earlier (its a vector that corresponds to which genes meet our filter criteria)
p7_vs_p5 <- p7_vs_p5[prefilter_p7_vs_p5,] #applying our filter we created earlier (its a vector that corresponds to which genes meet our filter criteria)
g7_vs_p7 <- g7_vs_p7[prefilter_g7_vs_p7,] #applying our filter we created earlier (its a vector that corresponds to which genes meet our filter criteria)

#Re-compute the library sizes for the normalization step next
g7_vs_g5$samples$lib.size <- colSums(g7_vs_g5$counts) #reseting the library size (lib.size) for the samples based on the number of genes being examined post-filtering
g7_vs_p5$samples$lib.size <- colSums(g7_vs_p5$counts)
p7_vs_p5$samples$lib.size <- colSums(p7_vs_p5$counts)
g7_vs_p7$samples$lib.size <- colSums(g7_vs_p7$counts)

#Calculate the normalization factors
g7_vs_g5 <- calcNormFactors(g7_vs_g5)
g7_vs_p5 <- calcNormFactors(g7_vs_p5)
p7_vs_p5 <- calcNormFactors(p7_vs_p5)
g7_vs_p7 <- calcNormFactors(g7_vs_p7)

#Look at the normalization factors
#This is an important step as we can see if one library is significantly bigger or smaller than another
g7_vs_g5$samples #normalization factors for glycerol pH 7 vs glycerol pH 5.7
g7_vs_p5$samples #normalization factors for glycerol pH 7 vs glycerol pH 5.7
p7_vs_p5$samples #normalization factors for pyruvate pH 7 vs pyruvate pH 5.7
g7_vs_p7$samples #normalization factors for glycerol pH 7 vs pyruvate pH 7

png('~/TranscriptomicsWorkshop/DiffExp/gly7vsgly5MDSplot.png')
plotMDS(g7_vs_g5)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr5MDSplot.png')
plotMDS(g7_vs_p5)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/pyr7vspyr5MDSplot.png')
plotMDS(p7_vs_p5)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr7MDSplot.png')
plotMDS(g7_vs_p7)
dev.off()

g7_vs_g5 <- estimateCommonDisp(g7_vs_g5, verbose=TRUE)
g7_vs_g5 <- estimateTagwiseDisp(g7_vs_g5)

g7_vs_p5 <- estimateCommonDisp(g7_vs_p5, verbose=TRUE)
g7_vs_p5 <- estimateTagwiseDisp(g7_vs_p5)

p7_vs_p5 <- estimateCommonDisp(p7_vs_p5, verbose=TRUE)
p7_vs_p5 <- estimateTagwiseDisp(p7_vs_p5)

g7_vs_p7 <- estimateCommonDisp(g7_vs_p7, verbose=TRUE)
g7_vs_p7 <- estimateTagwiseDisp(g7_vs_p7)

png('~/TranscriptomicsWorkshop/DiffExp/gly7vsgly5BCVplot.png')
plotBCV(g7_vs_g5, cex=1)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr5BCVplot.png')
plotBCV(g7_vs_p5, cex=1)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/pyr7vspyr5BCVplot.png')
plotBCV(p7_vs_p5, cex=1)
dev.off()
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr7BCVplot.png')
plotBCV(g7_vs_p7, cex=1)
dev.off()


et_gly7_vs_gly5 <- exactTest(g7_vs_g5)
et_gly7_vs_pyr5 <- exactTest(g7_vs_p5)
et_pyr7_vs_pyr5 <- exactTest(p7_vs_p5)
et_gly7_vs_pyr7 <- exactTest(g7_vs_p7)

top_gly7_vs_gly5 <- topTags(et_gly7_vs_gly5)
top_gly7_vs_pyr5 <- topTags(et_gly7_vs_pyr5)
top_pyr7_vs_pyr5 <- topTags(et_pyr7_vs_pyr5)
top_gly7_vs_pyr7 <- topTags(et_gly7_vs_pyr7)

summary(de_gly7_vs_gly5 <- decideTestsDGE(et_gly7_vs_gly5)) #show how many genes are significantly differentially expressed (downregulated = -1, upregulated = 1, no differential expression = 0)
detags <- rownames(g7_vs_g5)[as.logical(de_gly7_vs_gly5)] #create a vector of TRUE or FALSE using the gene names (rownames) on whether a gene is differentially expressed (TRUE) or not (FALSE)
png('~/TranscriptomicsWorkshop/DiffExp/gly7vsgly5scatter.png')
plotSmear(et_gly7_vs_gly5, de.tags=detags, cex=0.75) #plot the scatter plot
abline(h=c(-1, 1), col="blue") #the y-axis is log base 2 for fold-change so add a couple lines delineating 2-fold changes up and down
dev.off()

summary(de_gly7_vs_pyr5 <- decideTestsDGE(et_gly7_vs_pyr5))
detags <- rownames(g7_vs_p5)[as.logical(de_gly7_vs_pyr5)]
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr5scatter.png')
plotSmear(et_gly7_vs_pyr5, de.tags=detags, cex=0.75)
abline(h=c(-1, 1), col="blue")
dev.off()

summary(de_pyr7_vs_pyr5 <- decideTestsDGE(et_pyr7_vs_pyr5))
detags <- rownames(p7_vs_p5)[as.logical(de_pyr7_vs_pyr5)]
png('~/TranscriptomicsWorkshop/DiffExp/pyr7vspyr5scatter.png')
plotSmear(et_pyr7_vs_pyr5, de.tags=detags, cex=0.75)
abline(h=c(-1, 1), col="blue")
dev.off()

summary(de_gly7_vs_pyr7 <- decideTestsDGE(et_gly7_vs_pyr7))
detags <- rownames(g7_vs_p7)[as.logical(de_gly7_vs_pyr7)]
png('~/TranscriptomicsWorkshop/DiffExp/gly7vspyr7scatter.png')
plotSmear(et_gly7_vs_pyr7, de.tags=detags, cex=0.75)
abline(h=c(-1, 1), col="blue")
dev.off()

FDR_gly7_vs_gly5 <- p.adjust(et_gly7_vs_gly5$table$PValue, method="BH") #create a new vector of FDR corrected p-values
gly7_vs_gly5_outfile <- cbind(cpm(g7_vs_g5), et_gly7_vs_gly5$table, FDR_gly7_vs_gly5) #create a new dataframe with cpm (transcript abundances), summary data from the differential gene expression (exactTest()) table, and FDR corrected p-values
write.csv(gly7_vs_gly5_outfile, file="~/TranscriptomicsWorkshop/DiffExp/gly7_vs_gly5results.csv") #write the dataframe we just made out to a folder in Documents > example-RNAseq > analysis-results

FDR_gly7_vs_pyr5 <- p.adjust(et_gly7_vs_pyr5$table$PValue, method="BH")
gly7_vs_pyr5_outfile <- cbind(cpm(g7_vs_p5), et_gly7_vs_pyr5$table, FDR_gly7_vs_pyr5)
write.csv(gly7_vs_pyr5_outfile, file="~/TranscriptomicsWorkshop/DiffExp/gly7_vs_pyr5results.csv")

FDR_pyr7_vs_pyr5 <- p.adjust(et_pyr7_vs_pyr5$table$PValue, method="BH")
pyr7_vs_pyr5_outfile <- cbind(cpm(p7_vs_p5), et_pyr7_vs_pyr5$table, FDR_pyr7_vs_pyr5)
write.csv(pyr7_vs_pyr5_outfile, file="~/TranscriptomicsWorkshop/DiffExp/pyr7_vs_pyr5results.csv")

FDR_gly7_vs_pyr7 <- p.adjust(et_gly7_vs_pyr7$table$PValue, method="BH")
gly7_vs_pyr7_outfile <- cbind(cpm(g7_vs_p7), et_gly7_vs_pyr7$table, FDR_gly7_vs_pyr7)
write.csv(gly7_vs_pyr7_outfile, file="~/TranscriptomicsWorkshop/DiffExp/gly7_vs_pyr7results.csv")
