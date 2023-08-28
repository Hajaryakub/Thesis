#Code obtained from Liliane Conteville
#https://rpubs.com/lconteville/714853

library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")
library("igraph")

# Filter out singleton taxa from the phyloseq object
filtered_phylo_data<- filter_taxa(pruned, function (x) {sum(x > 0) > 1}, prune=TRUE) #down to 697

#using the object pruned(rarefied)  
pruned <- filter_taxa(filtered_phylo_data, function(x) {sum(x > 0) / length(x) >= 0.1}, prune=TRUE) #down to 383 

#transform the data to relative abundance and create a new phyloseq object:
relab_genera = transform_sample_counts(pruned, function(x) x / sum(x) * 100) 
head(otu_table(relab_genera)[,1:6])

# Calculate the Bray-Curtis distance among samples
bray_distance <- phyloseq::distance(relab_genera, method = "bray")

# Create the `bray_physeq` object directly
bray_physeq <- phyloseq(otu_table(relab_genera), sample_data(relab_genera), tax_table(relab_genera))

# Perform PCoA on the Bray-Curtis distance using the updated phyloseq object
pcoa <- ordinate(bray_physeq, method = "PCoA", distance= "bray")

# Plot the PCoA by group
plot_ordination(bray_physeq, pcoa, color="Group") + 
  geom_point(size=3) + 
  labs(color="Group") + 
  ggtitle("PCoA") + 
  theme_bw() 

# Plot the PCoA by Sex
plot_ordination(bray_physeq, pcoa, color = "Sex") +
  geom_point(size=3)+
  stat_ellipse(aes(group=Sex))

#plot the PCoA by Sex and Group
plot_ordination(bray_physeq, pcoa, color = "Group", shape="Sex") + 
  geom_point(size=4)+
  stat_ellipse(aes(group=Group))

#Perform PERMANOVA to see statisical difference
# Extract the necessary variables from the sample data
group_data <- sample_data(bray_physeq)$Group

# Create a data frame
sample_data_df <- data.frame(Group = group_data)

# Run adonis2 analysis
adonis2(bray_distance ~ Group, data = sample_data_df, permutations = 9999, method = "bray")

#Pr(>F) represents the p-value which is significant when < 0.05.
#not significant

#check the participant number to identify the outlier points
pcoa<- plot_ordination(bray_physeq, pcoa, color = "Group", shape="Sex") + 
  geom_point(size=4)+
  stat_ellipse(aes(group=Group))

# Extract the data for the points layer
point_data <- layer_data(pcoa, 1)

# Display the coordinates of each point
coordinates <- point_data[c("x", "y")]

# Filter the coordinates based on x-axis condition
outliers <- coordinates[coordinates$x < -0.3, ] #participant 19 is an outlier


#Removing outliers and replotting PCOA
#"BBK8038", "BBK8101", "BBK8109", "BBK8099"

outlier_samples <- c('BBK8038', 'BBK8101', 'BBK8109', 'BBK8099')

filtered_physeq <- bray_physeq
for (sample_id in outlier_samples) {
  filtered_physeq <- subset_samples(filtered_physeq, sample_names(filtered_physeq) != sample_id)
}

# Remove the samples from the relab_genera object
relab_genera_filtered <- subset_samples(relab_genera, !(sample_names(relab_genera) %in% outlier_samples))

#Perform PCoA on the updated phyloseq object
filtered_pcoa <- ordinate(filtered_physeq, method = "PCoA", distance = "bray")

#Replot the PCoA with the updated data
plot_ordination(filtered_physeq, filtered_pcoa, color="Group") + 
  geom_point(size=3) + 
  labs(color="Group") + 
  ggtitle("PCoA") + 
  theme_bw() 

# Plot the PCoA by Sex
plot_ordination(filtered_physeq, filtered_pcoa, color = "Sex") +
  geom_point(size=3)+
  stat_ellipse(aes(group=Sex))

# Create a new column "Age_category" based on age ranges
filtered_meta<-sample_data(filtered_physeq)
new_matrix<-as.matrix(filtered_meta)
filtered_meta<-as.data.frame(new_matrix)

age_metadata <- filtered_meta %>%
  mutate(Age_category = case_when(
    Age >= 18 & Age <= 25 ~ 1,
    Age >= 26 & Age <= 36 ~ 2,
    Age >= 37 & Age <= 47 ~ 3,
    TRUE ~ 4  # For ages higher than 47
  ))

#Update the metadata in your phyloseq object
sample_data(filtered_physeq)$Age_category <- age_metadata$Age_category

#Convert Age_category to a factor
sample_data(filtered_physeq)$Age_category <- factor(sample_data(filtered_physeq)$Age_category)

age_color_palette <- c("blue", "red", "orange", "purple")  # Add more colors if needed

# Use scale_color_manual to set colors based on Age_category values
plot_ordination(filtered_physeq, filtered_pcoa, color = "Age_category") +
  geom_point(size = 3) +
  scale_color_manual(values = age_color_palette)+
  stat_ellipse(aes(group=Age_category))


#Re-Run Permanova to see statisical difference
# Re-Calculate the Bray-Curtis distance among samples
bray_distance_filtered <- phyloseq::distance(relab_genera_filtered, method = "bray")

#Extract the necessary variables from the sample data
filtered_group_data <- sample_data(filtered_physeq)$Group

# Create a data frame
sample_data_df2 <- data.frame(Group = filtered_group_data)

# Run adonis2 analysis
adonis2(bray_distance_filtered ~ Group, data = sample_data_df2, permutations = 9999, method = "bray") ##Still not significant##


#Extract the necessary variables from the sample data to test significance of Sex
filtered_sex_data <- sample_data(filtered_physeq)$Sex

# Create a data frame
sample_data_sex <- data.frame(Sex = filtered_sex_data)

# Run adonis2 analysis
adonis2(bray_distance_filtered ~ Sex, data = sample_data_sex, permutations = 9999, method = "bray")
#significant 1e-04

#Extract the necessary variables from the sample data to test significance of Age
filtered_age_data <- sample_data(filtered_physeq)$Age_category

# Create a data frame
sample_data_age <- data.frame(Age = filtered_age_data)

# Run adonis2 analysis
adonis2(bray_distance_filtered ~ Age, data = sample_data_age, permutations = 9999, method = "bray")
#1e-04

#heatmaps
nooutlier_physeq <- subset_samples(pruned, !rownames(sample_data(pruned)) %in% c('BBK8038', 'BBK8101', 'BBK8109', 'BBK8099'))
TopNGenus <- names(sort(taxa_sums(nooutlier_physeq), TRUE)[1:4]) 
Top4Genus <- prune_taxa(TopNGenus, nooutlier_physeq)
plot_heatmap(Top4Genus)
plot_heatmap(Top4Genus, "PCoA", "bray", low="#66CCFF", high="#000033")
