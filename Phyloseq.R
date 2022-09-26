#########PACKAGES############
library("dplyr")
library("phyloseq")
library("ggplot2")
library("ggpubr")
library("microbiome")
library("gganimate")
library("viridis")
library("rstatix")

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
##### SUBSETTING DATA #####
############################

#If you want to specifically subset your sample for the REMAINDER of the analysis, you can do that by patient or whatever you choose
#phylo_obj_f<- subset_samples(phylo_obj_f, Type=="Abscess")


############################
##### Filter Singletons #####
############################
phylo_obj_f_pruned_fs <- prune_taxa(taxa_sums(phylo_obj_f_pruned) > 1, phylo_obj_f_pruned) 
summarize_phyloseq(phylo_obj_f_pruned_fs)

#Use singleton Filtered Data for all other analysis!
#Write out a file of the fs data so you can look at it in excel or whatever you want to do 
phylo_obj_f_pruned_fs_df <- psmelt(phylo_obj_f_pruned_fs)
write.csv(phylo_obj_f_pruned_fs_df,"/Users/kriegema/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/FusoZinProt/Counts-Fs.csv", row.names = FALSE)


############################
#### PREPROCESSING DATA ####
############################
#translate your data by adding pseudo count
phylo_obj_f_pruned_fs_t = transform_sample_counts(phylo_obj_f_pruned_fs, function(x) (x + 0.000001 - min(x))) 
otu_table(phylo_obj_f_pruned_fs_t)

#transform to relative abundance with translated pseudo count
phylo_obj_f_pruned_fs_t_n = transform_sample_counts(phylo_obj_f_pruned_fs_t, function(x) (x)/sum(x)) ##can use microbiome package instead of this, see line below
otu_table(phylo_obj_f_pruned_fs_t_n)

#microbiome::transform(phylo_obj_fs_t, transform=clr) ##Can also do this filtering step insetad of the transform_sample_counts step above. Kris says that he prefers CLR, but the log10 transformation is widely accepted so we need to do this for transformations to publish

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
#Species
phylo_obj_f_pruned_fs_t_n_sglom_df$Species <- gsub("s__","",as.character(phylo_obj_f_pruned_fs_t_n_sglom_df$Species))

