library(ALDEx2)
packageVersion("ALDEx2") #‘1.30.0’

pruned<-nooutlier_physeq

# Prepare your `phyloseq` object ("pruned") and subset it to only two conditions (post and pre supplement) for ALDEx2 to run
subset_physeq <- subset_samples(pruned, Group %in% c("Pre_ERME", "Post_ERME"))
# Extract the count data from phyloseq object
count_data <- otu_table(subset_physeq)

# Convert the count data to a matrix
count_matrix <- as.matrix(count_data)

# Extract the group or condition variable from sample data
group_var <- sample_data(subset_physeq)$Group

# Perform differential abundance analysis with ALDEx2
aldex_data<-aldex(count_matrix, conditions = group_var, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE)

# Order results by effect size
ordered_results <- aldex_data[order(aldex_data$effect, decreasing = TRUE), ]

#MA and Effect plots of ALDEx2 output
#In both plots features that are not significant are in grey or black. 
par(mfrow=c(1,2))
aldex.plot(aldex_data, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(aldex_data, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

##new subset to compare post and pre placebo##
subset_placebo <- subset_samples(pruned, Group %in% c("Pre_Placebo", "Post_Placebo"))
# Extract the count data from phyloseq object
count_data2 <- otu_table(subset_placebo)

# Convert the count data to a matrix
count_matrix2 <- as.matrix(count_data2)

# Extract the group or condition variable from sample data
group_var2 <- sample_data(subset_placebo)$Group

# Perform differential abundance analysis with ALDEx2
aldex_data2<-aldex(count_matrix2, conditions = group_var2, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE)

#MA and Effect plots of ALDEx2 output
#In both plots features that are not significant are in grey or black. 
par(mfrow=c(1,2))
aldex.plot(aldex_data2, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(aldex_data2, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

#BOTH GROUP COMPARISONS HAVE NO SIGNIFICANCE DIFFERENCE

#plotting a network on subset of only PRE and POST ERME
plot_net(subset_physeq, distance = "bray", type = "samples", maxdist = 0.7, color = "Group", shape = "Group")
ig <- make_network(subset_physeq, max.dist=0.8)
plot_network(ig, subset_physeq, color="Group", shape="Group")


#Save files for Lefse Analysis
pruned_otus<- otu_table(pruned)
pruned_metadata<-sample_data(pruned)
pruned_taxons<- tax_table(pruned)

write.csv(pruned_otus, file = "/rds/homes/h/hyb288/Thesis Data/otu_table_lefse.csv", row.names = TRUE)
write.csv(pruned_taxons, file = "/rds/homes/h/hyb288/Thesis Data/taxonomy_lefse.csv", row.names = TRUE)
write.csv(pruned_metadata, file = "/rds/homes/h/hyb288/Thesis Data/metadata_lefse.csv", row.names = TRUE)

#Check DAA based on Gender and group: pre erme and post erme
sex_var <- sample_data(subset_physeq)$Sex

# Perform differential abundance analysis with ALDEx2
aldex_sex<-aldex(count_matrix, conditions = sex_var, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE)

# Order results by effect size
ordered_results_sex <- aldex_sex[order(aldex_sex$effect, decreasing = TRUE), ]

#MA and Effect plots of ALDEx2 output
#In both plots features that are not significant are in grey or black. 
par(mfrow=c(1,2))
aldex.plot(aldex_sex, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(aldex_sex, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

#More plots showing median log2 abundance. Change sex_var to group_var
x <- aldex.clr(count_matrix, group_var, mc.samples=16, denom="all", verbose=F)
x.tt <- aldex.ttest(x, paired.test=FALSE, verbose=FALSE)
x.effect <- aldex.effect(x, CI=T, verbose=FALSE, paired.test=FALSE)
x.all <- data.frame(x.tt,x.effect)
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")


library(ANCOMBC)
packageVersion("ANCOMBC")
library(DT)
set.seed(123)
out = ancombc(data = NULL, assay_name = NULL,
              tax_level = "family", phyloseq = pruned,
              formula = "Age + Sex + Bodyweight.kg + Group",
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
              group = "Group", struc_zero = TRUE, neg_lb = FALSE,
              tol = 1e-5, max_iter = 100, conserve = TRUE,
              alpha = 0.05, global = TRUE, n_cl = 1, verbose = TRUE)
res = out$res
res_global = out$res_global

col_name = c("Taxon", "Intercept", "Age", "Sex", "Bodyweight.kg", 
             "Pre_ERME-Post_ERME", "")

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result") 

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)