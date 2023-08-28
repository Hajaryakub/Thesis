rm(list=ls())
options(bitmapType = 'cairo')
install.packages("dplyr")  # Install the latest version of dplyr

library(phyloseq) ; packageVersion("phyloseq") #1.42.0
library(vegan) ; packageVersion("vegan") #2.6.2
library(tidyverse) ; packageVersion("tidyverse") #1.3.1
library(dplyr) ; packageVersion("dplyr") #1.0.8
library(ggplot2) ; packageVersion("ggplot2") #3.3.5
library(DESeq2) ; packageVersion("DESeq2") #1.38.3
library(dendextend) ; packageVersion("dendextend") #1.15.2

setwd("/rds/homes/h/hyb288/Thesis Data")

#read the table containing species, taxonomy, and sample counts
profile<- read.table("profiles.tsv", header=TRUE, sep = "\t")

#create a separate taxonomy table from the profiles
taxonomy<- profile %>%
  separate(GTDB.taxonomy, into =c("domain", "phylum", "class", "order", "family", "genus"), sep = ";", extra = "drop", fill = "right") %>%
  dplyr::select(-c(8:76)) 

# Change the column name from "Taxon" to "species" using base R indexing
colnames(taxonomy)[colnames(taxonomy) == "Taxon"] <- "species"

#removing the "d_, p_" etc from each class column
taxonomy$species<-gsub("s__","",as.character(taxonomy$species))
taxonomy$domain<-gsub("d__","",as.character(taxonomy$domain))
taxonomy$phylum<-gsub("p__","",as.character(taxonomy$phylum))
taxonomy$class<-gsub("c__","",as.character(taxonomy$class))
taxonomy$order<-gsub("o__","",as.character(taxonomy$order))
taxonomy$family<-gsub("f__","",as.character(taxonomy$family))
taxonomy$genus<-gsub("g__","",as.character(taxonomy$genus))

#creating the otu table by deleting unneccessary columns
otu_data<- dplyr::select(profile, Taxon, everything())
otu_data<- dplyr::select(otu_data, -GTDB.taxonomy) 
# Change the column name from "Taxon" to "species" using base R indexing
colnames(otu_data)[colnames(otu_data) == "Taxon"] <- "species"
otu_data$species<-gsub("s__","",as.character(otu_data$species))

#set the rownames
rownames(otu_data)<- otu_data$species
rownames(taxonomy)<-taxonomy$species
# Remove the extra "species" column from otu data, and taxonomy
otu_data$species <- NULL
taxonomy$species<- NULL

#read the metadata
metadata<- read.table("MetaData.txt", header = TRUE, sep = "\t", row.names=1)

# Reorder the metadata based on the sample order in the otu table
sample_order<- match(colnames(otu_data), rownames(metadata))
metadata<- metadata[sample_order, ]

# Reorder the taxonomy data based on the species order in the otu table
species_order<- match(rownames(otu_data), rownames(taxonomy))
taxonomy<- taxonomy[species_order, ]

#save the output files
write.csv(otu_data, file = "/rds/homes/h/hyb288/Thesis Data/otu_table_fromR.csv", row.names = TRUE)
write.csv(taxonomy, file = "/rds/homes/h/hyb288/Thesis Data/taxonomy_fromR.csv", row.names = TRUE)
write.csv(metadata, file = "/rds/homes/h/hyb288/Thesis Data/metadata_fromR.csv", row.names = TRUE)

#creating the phyloseq object
taxonomy <- as.matrix(taxonomy)
physeq <- phyloseq(otu_table(otu_data, taxa_are_rows = TRUE), sample_data(metadata), tax_table(taxonomy))

#checking the sequencing depth
sequencing_depth <- sample_sums(physeq)

library(ggplot2)
# Create a data frame with sample names and sequencing depth
depth_df <- data.frame(Sample = sample_names(physeq), Depth = sequencing_depth)

# Plot the sequencing depth using ggplot2
ggplot(depth_df, aes(x = Sample, y = Depth)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Sample") +
  ylab("Sequencing Depth") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#distribution of taxa
taxa_sums(physeq)

#checking singletons
otu_occurrences <- rowSums(otu_data > 0)
singletons <- names(otu_occurrences[otu_occurrences == 1]) #173 singletons

#Checking pcoa (bray-curtis) before normalization
# Create the `bray_physeq` object directly
brays_data <- phyloseq(otu_table(physeq), sample_data(physeq), tax_table(physeq))

# Perform PCoA on the Bray-Curtis distance using the updated phyloseq object
first_pcoa <- ordinate(brays_data, method = "PCoA", distance= "bray")

plot_ordination(brays_data, first_pcoa, color="Group") + 
  geom_point(size=1) + labs(col="Group") + 
  geom_text(aes(label=rownames(metadata), hjust=0.3, vjust=-0.4)) + 
  ggtitle("PCoA") + 
  theme_bw() + theme(legend.position="none")


###################################################
#Plotting rarefaction curve
#run ggrare_function.R; output is ps.rarefied
###################################################


#Abundance plot at different taxonomic levels
#https://github.com/joey711/phyloseq/issues/864
#plot_bar(ps.rarefied, fill="family") + facet_wrap(~Group, scales="free_x", nrow=1)
library(RColorBrewer)
taxonomy_df<- as.data.frame(taxonomy)

Phyla_num <- length(levels(as.factor(taxonomy_df$phylum)))
getPalette<- colorRampPalette(brewer.pal(9, "Set2")) 
PhylaPalette<- getPalette(Phyla_num)
pd <- psmelt(vst_physeq)

abundance_plot<- ggplot(pd, aes(x = Sample, y = Abundance, factor(phyla), fill = factor(phylum))) + 
  geom_bar(stat = "identity") + facet_wrap(~Group, scales = "free_x") + scale_fill_manual(values = PhylaPalette) + 
  guides(fill=guide_legend(ncol=2)) 

ggsave("abundance_plot_phyla.png", abundance_plot, width = 30, height = 12)

###################################################
#run Alpha_diversity.R 
###################################################

###################################################
#run Beta_diversity.R 
###################################################

###################################################
#Aldex2 and ANCOM-BC
#run DAA
###################################################
