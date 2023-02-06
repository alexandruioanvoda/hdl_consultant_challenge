# Downloading and reading in the data -------
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz",
              "GSE40279_series_matrix.txt.gz")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/suppl/GSE40279_average_beta.txt.gz",
              "GSE40279_average_beta.txt.gz")

library(data.table)
bmat <- as.data.frame(fread("GSE40279_average_beta.txt.gz"))
#View(head(bmat))
x <- readLines(gzfile("GSE40279_series_matrix.txt.gz"))
#View(as.data.frame(x))

# Boilerplate to obtain the metadata matrix -------
meta <- data.frame(sample = strsplit(x = x[which(startsWith(x, "!Sample_source_name_ch1"))[1]], split = "\t"),
                   age = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[1]], split = "\t"),
                   source = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[2]], split = "\t"),
                   plate = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[3]], split = "\t"),
                   gender = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[4]], split = "\t"),
                   ethnicity = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[5]], split = "\t"),
                   tissue = strsplit(x = x[which(startsWith(x, "!Sample_characteristics_ch1"))[6]], split = "\t"))
# Ok, let's clean the meta
meta[1:5,]
meta <- meta[-1,]
colnames(meta) <- c("sample", "age", "source", "plate", "gender", "ethnicity", "tissue")
meta$sample <- gsub('"', '', meta$sample, fixed=TRUE)
meta$age <- gsub('age \\(y\\)\\: |"', '', meta$age)
meta$source <- gsub('source\\: |"', '', meta$source)
meta$plate <- gsub('plate\\: |"', '', meta$plate)
meta$gender <- gsub('gender\\: |"', '', meta$gender)
meta$ethnicity <- gsub('ethnicity\\: |"', '', meta$ethnicity)
meta$tissue <- gsub('tissue\\: |"', '', meta$tissue)
meta[1:5,]
meta$age <- as.numeric(meta$age)
meta$plate <- as.numeric(meta$plate)

hist(meta$age) # minimum age is 19, max is 101, seems like a normal distribution with mean of 60 or so, and a standard deviation of 10-20 years or so

table(meta$gender) # roughly 50-50, that's good
table(meta$ethnicity) # roughly two thirds euro/caucasian and one third hispanic/mexican
table(meta$tissue) # so, just one tissue (whole blood), that's good.
table(meta$plate) # plate 10 seems to be a small outlier, with just 17 samples (and plate 9 with 48 samples), but the others have 80-90 samples
table(meta$source) # Boston provides the least samples (35)

all(meta$sample %in% colnames(bmat)) # TRUE # Good.

# let's clean up intermediary variables
rm(x); gc()

# Boilerplate to reformat the beta matrix -------
## Ok, now I need to reformat bmat (the Beta Matrix) into something that has variables (CpG island IDs) as columns and samples (individuals) as rows
#View(head(bmat))
rownames(bmat) <- bmat[,1]
bmat <- bmat[,-1]
bmat <- t(bmat)
dim(bmat) # So 656 samples (individuals) and 473034 variables (probes)
class(bmat)

### Now, let's add the metadata to this.
df <- cbind(bmat, meta[match(rownames(bmat), meta$sample),])
# saveRDS(df, "df_GSE40279.rds") # Saving a checkpoint of all the data here.
# df <- readRDS("df_GSE40279.rds")
#View(df[1:10,c("cg09809672", "age")])
plot(cg09809672~age, data = df) # Specifically looking at this one because I selected it from https://www.frontiersin.org/articles/10.3389/fgene.2021.759357/full
#View(df[1:100, 1:100])

# Splitting the data into two sets: for training and for testing -------
set.seed(42) # No, I have not read the Hitchhiker's Guide to the Galaxy. This is just my seed number.
train_idx <- sample(1:nrow(df), size = as.integer(0.75*nrow(df)))

train_df <- df[train_idx,]
test_df <- df[-train_idx,]

dim(train_df)
dim(test_df)
# Let's clean up intermediary variables
rm(train_idx, df); gc()
#rm(bmat,meta); gc()

# colnames(df)[473035:473041-7] # So as expected columns before are CpG islands
# colnames(df)[473035:473041] # And columns after are sample metadata

# Actually building the elastic net model and saving the weights -------
library(glmnet)
library(dplyr)
set.seed(42)

system.time({model <- glmnet(train_df %>% select(starts_with("cg")),
                             train_df$age,
                             alpha=0.1,
                             lambda=1)}) # took 24 seconds
saveRDS(model, "model_task1.rds")

# Evaluating the model -------
# Short check that one important CpG has a weight (aka it isn't squashed by the regularization)
model[["beta"]]["cg09809672",]

# Compute predicted values
predicted_ages <- predict(model, as.matrix(test_df %>% select(starts_with("cg"))))

residuals <- predicted_ages - test_df$age # Compute the residuals
total_variance <- var(test_df$age) * (length(test_df$age) - 1) # Compute the total variance
residual_variance <- sum(residuals^2) # Compute the residual variance
r_squared <- 1 - (residual_variance / total_variance) # Compute R-squared
mae <- mean(abs(predicted_ages - test_df$age)) # Compute MAE
#rmse <- sqrt(mean((residuals)^2)) # Compute RMSE

print(c(r_squared, mae))

plot(test_df$age, predicted_ages, xlab="Actual Age", ylab="Predicted Age")
#abline(0, 1, col="red")

