#########PACKAGES############
library("dplyr")
library("phyloseq")
library("ggplot2")

########################
#### IMPORTING DATA ####
########################
getwd()
#Importing the .biom file from Qiime2 output
#The ASV biom file is generated in excel and converted back to .biom format using Qiime2
ps1 = "Qiime2/Trimmed/Files_For_Phyloseq/feature_table_w_taxonomy.biom" #The biom should not have been filtered for singletons! this will skew downstream analysis
myData_ps1 = import_biom(ps1) 
myData_ps1

metadata = import_qiime_sample_data("FusoZinProt/mapfile.txt")
head(metadata)

RootedTree = "FusoZinProt/Qiime2/Trimmed/Tree/Unfiltered_Rooted_tree_for_phyloseq/tree.nwk"
Tree <- read_tree(RootedTree)
Tree

#Add the Taxonomic Ranks to the different Taxonomic levels
myData_ps1
colnames(tax_table(myData_ps1)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rank_names(myData_ps1)

#Generate the Phyloseq Object
phylo_obj = merge_phyloseq(myData_ps1,metadata,Tree)
phylo_obj

check <- data.frame(sample_data(phylo_obj))
check


#########################
#### FILTERING DATA ####
#########################

#Filter Unassigned
phylo_obj_f <- subset_taxa(phylo_obj, Phylum != "Unassigned")
summarize_phyloseq(phylo_obj_f)


##############
###READ COUNTS
##############

sample_sums(phylo_obj)
sort(sample_sums(phylo_obj))
hist(sample_sums(phylo_obj), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=12)

metadata$total_reads <- sample_sums(phylo_obj)
read_counts <- arrange(metadata, total_reads) ##COULD EXPORT THIS FOR INFO

#Remove values with less than 10,000 total reads:
phylo_obj_f_pruned <- prune_samples(sample_sums(phylo_obj_f)>=10000, phylo_obj_f)
metadata_pruned <- data.frame(sample_data(phylo_obj_f_pruned))
metadata_pruned$total_reads <- sample_sums(phylo_obj_f_pruned)
read_counts_pruned <- arrange(metadata_pruned, total_reads)
summarize_phyloseq(phylo_obj_f_pruned)

############################
##### Filter Singletons #####
############################
phylo_obj_f_pruned_fs <- prune_taxa(taxa_sums(phylo_obj_f_pruned) > 1, phylo_obj_f_pruned) 
summarize_phyloseq(phylo_obj_f_pruned_fs)

#Use singleton Filtered Data for all other analysis!
#Write out a file of the fs data so you can look at it in excel or whatever you want to do 
phylo_obj_f_pruned_fs_df <- psmelt(phylo_obj_f_pruned_fs)
write.csv(phylo_obj_f_pruned_fs_df,"FusoZinProt/Counts-Fs.csv", row.names = FALSE)


############################
#### PREPROCESSING DATA ####
############################
#translate your data by adding pseudo count
phylo_obj_f_pruned_fs_t = transform_sample_counts(phylo_obj_f_pruned_fs, function(x) (x + 0.000001 - min(x))) 
otu_table(phylo_obj_f_pruned_fs_t)

#transform to relative abundance with translated pseudo count
phylo_obj_f_pruned_fs_t_n = transform_sample_counts(phylo_obj_f_pruned_fs_t, function(x) (x)/sum(x))
otu_table(phylo_obj_f_pruned_fs_t_n)


#########################
### Agglomerate Data ####
#########################
phylo_obj_f_pruned_fs_t_n_sglom <- tax_glom(phylo_obj_f_pruned_fs_t_n, taxrank = "Species", NArm = TRUE) #Species level


#########################
#### CONVERTING DATA ####
#########################
#convert to data.frame using phyloseq psmelt function
phylo_obj_f_pruned_fs_t_n_sglom_df <- psmelt(phylo_obj_f_pruned_fs_t_n_sglom)

#Clean dataframe so that we don't have the s__ in front of Species and so there are no "_" characters
phylo_obj_f_pruned_fs_t_n_sglom_df$Species <- gsub("s__","",as.character(phylo_obj_f_pruned_fs_t_n_sglom_df$Species))

                                                    
################
#ggplots Bar Plot
################
#Make sure all the abundances are 1
aggregate(phylo_obj_f_pruned_fs_t_n_sglom_df$Abundance, list(phylo_obj_f_pruned_fs_t_n_sglom_df$Sample), FUN=sum)
phylo_obj_f_pruned_fs_t_n_sglom_df


phylo_obj_f_pruned_fs_t_n_sglom_df[phylo_obj_f_pruned_fs_t_n_sglom_df == "nucleatum subsp. animalis"] <- "animalis"
phylo_obj_f_pruned_fs_t_n_sglom_df[phylo_obj_f_pruned_fs_t_n_sglom_df == "nucleatum subsp. nucleatum"] <- "nucleatum"
phylo_obj_f_pruned_fs_t_n_sglom_df[phylo_obj_f_pruned_fs_t_n_sglom_df == "nucleatum subsp. fusiforme"] <- "vincentii"
phylo_obj_f_pruned_fs_t_n_sglom_df[phylo_obj_f_pruned_fs_t_n_sglom_df == "nucleatum subsp. polymorphum"] <- "polymorphum"
phylo_obj_f_pruned_fs_t_n_sglom_df[phylo_obj_f_pruned_fs_t_n_sglom_df == "nucleatum periodonticum"] <- "periodonticum"

#For everything lumped together
phylo_obj_f_pruned_fs_t_n_sglom_df%>%
  ggplot(aes(x = Sample, y = Abundance, fill = Species, reorder_within(Species=="nucleatum subsp. animalis", Abundance)))+
  scale_x_discrete(drop = TRUE)+
  geom_bar(aes(), stat="identity", size=1)+
  theme(strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"))+
  theme(axis.text.x=element_blank()) +
  theme(legend.text=element_text(face="italic"))+ 
  theme(strip.placement="top")+
  facet_grid(~factor(Type, levels=c("Plaque", "Abscess")), scales="free") +
  scale_fill_discrete(limits = c("animalis", "nucleatum", "perdiodonticum", "polymorphum", "vincentii"))+
  scale_fill_manual(values=c("#5290f2", "#F9CBAD", "#969595", "#f5f573", "#4ae059"))

ggsave(path = "~/Desktop/", filename="RelAbund.pdf", height =3.5, device='pdf', dpi=1000)
                                                    
                                                    
################
#CORRELATION MATRIX
################
library("ggpubr")
library("Hmisc")
library(tidyverse)
library(corrplot)

Counts_Fs_relabundance_paired_reformatted_forCor = subset(Counts_Fs_relabundance_paired_reformatted_forR, select = -c(1)) #remove the first row from the original df
Counts_Fs_relabundance_paired_reformatted_forCor <- Counts_Fs_relabundance_paired_reformatted_forCor %>% remove_rownames %>% column_to_rownames(var="Patient") #turn the first row (patient) into row names

Counts_Fs_relabundance_paired_reformatted_forCor2 <-
  Counts_Fs_relabundance_paired_reformatted_forCor %>%
  rename_all(funs(c("FFAbscess_animalis", "JJAbscess_fusiforme", "GGAbscess_nucleatum", "IIAbscess_polymorphum", "HHAbscess_periodonticum",
                    "AAPlaque_animalis", "EEPlaque_fusiforme", "BBPlaque_nucleatum", "DDPlaque_polymorphum", "CCPlaque_periodonticum"))) %>%
  select(AAPlaque_animalis, BBPlaque_nucleatum, CCPlaque_periodonticum, DDPlaque_polymorphum, EEPlaque_fusiforme, 
         FFAbscess_animalis, GGAbscess_nucleatum, HHAbscess_periodonticum, IIAbscess_polymorphum, JJAbscess_fusiforme)

Counts_Fs_relabundance_paired_reformatted_forCor2

Res = rcorr(as.matrix(Counts_Fs_relabundance_paired_reformatted_forCor2), type="spearman")
Res$P

res2 = cor(Counts_Fs_relabundance_paired_reformatted_forCor2)
testRes2 = cor.mtest(Counts_Fs_relabundance_paired_reformatted_forCor2,conf.level = 0.95)
testRes2

plot.new()
corrplot(res2, p.mat = testRes2$p, sig.level = 0.05, 
         type = "upper", insig='blank', addCoef.col ='black', 
         order = "hclust", tl.cex = .6,
         number.cex = 0.6, tl.col = "black", tl.srt = 45, col = COL2('BrBG'))


png(width=7, height=5, units="in", res=800, file="~/Desktop/Cor-matrix.png")
corrplot(res2, p.mat = testRes2$p, 
         insig = 'label_sig', 
         sig.level = c(0.001, 0.01, 0.05), 
         type = "upper", 
         tl.cex = .6, 
         order = "alphabet",
         pch.cex = 0.9, 
         pch.col = 'grey20',
         tl.col = "black", 
         tl.srt = 45)
dev.off()                                                    
                                                    
                                                  

