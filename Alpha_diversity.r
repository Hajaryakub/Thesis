#Code obtained from Liliane Conteville
#https://rpubs.com/lconteville/714853


#Calculating Alpha diversity
#trimming OTUS not present in any samples
pruned<- prune_species(speciesSums(ps.rarefied) >0, ps.rarefied)
pruned = transform_sample_counts(pruned, round)

#calculate richness values
richness <- estimate_richness(pruned)
head(richness)

#check if the shannon diversity is normally distributed
hist(richness$Shannon, main="Shannon index", xlab="") # it is not

#For Non-normally distributed data, we can use Kruskal-Wallis Rank Sum Test to check if the variable Group impacts the Shannon diversity
kruskal.test(richness$Shannon ~ sample_data(pruned)$Group) #no significant difference

#get a list with the p-values resulted of the Wilcoxon Tests considering each pair of groups:
pairwise.wilcox.test(richness$Shannon, sample_data(pruned)$Group, p.adj = "bonf")

#check if the observed diversity is normally distributed
hist(richness$Observed, main="Observed index", xlab="") # it is

#Since it is, we can run anova tests and check if the variable Group impacts the Shannon diversity:
anova.ob = aov(richness$Observed ~ sample_data(pruned)$Group)
summary(anova.ob)

#Based on the anova tests results, it’s possible to compute the Tukey Honest Significant Differences:
TukeyHSD(anova.ob)

#check if the Chao1 diversity is normally distributed
hist(richness$Chao1, main="Chao1 index", xlab="") #it is

#Since it is, we can run anova tests and check if the variable Group impacts the Shannon diversity:
anova.ch = aov(richness$Chao1 ~ sample_data(pruned)$Group)
summary(anova.ch)

#Based on the anova tests results, it’s possible to compute the Tukey Honest Significant Differences:
TukeyHSD(anova.ch)


#plotting all alpha diversity measures
#plot_richness(pruned, x="Group") + geom_boxplot()

#plotting shannon diversity 
plot_richness(pruned, x="Group", measures="Shannon", color = "Group")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16)) +
  ggtitle("Shannon") +
  theme(plot.title = element_text(size = 20))

#plotting observed diversity
plot_richness(pruned, x="Group", measures="Observed", color = "Group")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16)) +
  ggtitle("Observed") +
  theme(plot.title = element_text(size = 20))

#plotting Chao1 diversity
plot_richness(pruned, x="Group", measures="Chao1", color = "Group")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=16))+
  ggtitle("Chao1") +
  theme(plot.title = element_text(size = 20))
#no significant differences in alpha diversity between any of the groups