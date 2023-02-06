download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/suppl/GSE42861_processed_methylation_matrix.txt.gz",
              "GSE42861_processed_methylation_matrix.txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz",
              "GSE42861_series_matrix.txt.gz")

library(data.table)
bmat <- as.data.frame(fread("GSE42861_processed_methylation_matrix.txt.gz"))
#View(head(bmat))
x <- readLines(gzfile("GSE42861_series_matrix.txt.gz"))
#View(as.data.frame(x))

# Boilerplate to obtain the metadata matrix -------
meta <- data.frame(sample = strsplit(x = x[which(startsWith(x, "!Sample_supplementary_file"))[1]], split = "\t"),
                   age = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[4]], split = "\t"),
                   gender = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[5]], split = "\t"),
                   smoker = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[6]], split = "\t"))
# Ok, let's clean the meta
meta[1:5,]
meta <- meta[-1,]
colnames(meta) <- c("sample", "age", "gender", "smoker")
meta$sample <- sapply(meta$sample, function(i) {
  paste(strsplit(i, "_")[[1]][2:3], collapse = "_")
})
meta$age <- gsub('age\\: |"', '', meta$age)
meta$gender <- gsub('gender\\: |"', '', meta$gender)
meta$smoker <- gsub('smoking status\\: |"', '', meta$smoker)
meta[1:5,]
meta$age <- as.numeric(meta$age)

hist(meta$age) # minimum age is 19, max is 101, seems like a normal distribution with mean of 60 or so, and a standard deviation of 10-20 years or so
table(meta$gender) # at least twice more women than men in the dataset
table(meta$smoker) # there's 2 NAs in the data
all(meta$sample %in% colnames(bmat)) # TRUE # Good.

# let's clean up intermediary variables
rm(x); gc()

# Boilerplate to reformat the beta matrix -------
## Ok, now I need to reformat bmat (the Beta Matrix) into something that has variables (CpG island IDs) as columns and samples (individuals) as rows
#View(head(bmat))
rownames(bmat) <- bmat[,1]
bmat <- bmat[,-1]
bmat <- t(bmat)
dim(bmat) # So 690 samples (individuals) and 485577 variables (probes)
class(bmat)

# Now, let's add the metadata to this.
df <- cbind(bmat, meta[match(rownames(bmat), meta$sample),])
df$gender <- factor(df$gender)
plot(cg20926353~gender, data = df) # Taken from https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-014-0040-6
#View(df[1:100, 1:100])

ncol_bmat <- ncol(bmat)
rm(bmat, meta); gc()

sumstats <- as.data.frame(do.call("rbind", pbapply::pblapply(1:ncol_bmat, function(i) {
  return(c(summary(lm(df[,i] ~ df$gender))[["coefficients"]][2,],
           colnames(df)[i]))
}))) # takes about 2 hours.

for (i in 1:4) {sumstats[,i] <- as.numeric(sumstats[,i])}
sumstats <- sumstats[,-3]
colnames(sumstats) <- c("Effect", "StdErr", "P", "Probe")
sumstats$Padj <- p.adjust(sumstats$P, method="bonferroni")
table(sumstats$Padj < 0.05)
table(sumstats$P < 0.05)
saveRDS(sumstats, "sumstats_task2.rds")

p <- ggplot(sumstats, aes(x=Effect, y=-log10(P), color = Padj < 0.05)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05/nrow(sumstats)), color="blue", linetype="dashed", alpha=0.5)
ggsave("plot.png", plot=p, width=1700, height=1300, units="px")

library(methylGSA)

cpg.pval <- sumstats$P
names(cpg.pval) <- sumstats$Probe
cpg.pval[cpg.pval==0] <- 9.881313e-324 # min(cpg.pval[cpg.pval!=0]) # Because methylGSA does not accept P=0.
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
resgo = methylglm(cpg.pval = cpg.pval, minsize = 200, 
                 maxsize = 500, GS.type = "GO")
resreact = methylglm(cpg.pval = cpg.pval, minsize = 200, 
                  maxsize = 500, GS.type = "Reactome")
barplot(resreact, num = 10, colorby = "pvalue")



