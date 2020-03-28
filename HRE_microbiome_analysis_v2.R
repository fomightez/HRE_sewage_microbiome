##=======================================================================##
##=======================================================================##
#=========== Microbiome analysis of Piermont Pier samples 2017-18 ========#
##=======================================================================##
#================= From Hudson River Estuary ("HRE") =====================#
#=== Samples submitted for sequencing to MSU Genomics Core in July 2019 ==#
# ASVs were identified by qiime2 (See my notebook- "HRE_qiime2_notes.md") #
##=======================================================================##
##=======================================================================##
##=======================================================================##
# This analysis started July 22, 2019

# Feb 2020- decided to make "v2" of this script based on new needs/ discussion with Elise
# Switching from using vegan to using packages from phyloseq 
# because phyloseq allows for me to do the differential abundance analysis
# Also there are some nice workflows that I can adapt from the internet


# Generally began by folowing this workflow
# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html

# Some installations- this got tricky

# install.packages("knitr")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")

# NOTE- had to play around with this because I am using a newer R version. 
# Found some helpful info here https://bioconductor.org/install/


## Load packages

library(stats)
library(tidyverse)
library(vegan)
library(ggthemes)
library(reshape2)
library(devtools)
library(viridis)
library("knitr")
library("BiocManager")

.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  install(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)





# rm(list=ls(all=TRUE))



# Skip ahead in tutorial to "Combine data in phyloseq object"
# because I already did the analysis in qiime and have a count table, etc
# I can just import those instead of starting from scratch

# Note- I copied some output files from qiime in a directory called "SeqAnalysis_mod"

# Found this helpful pdf online: https://github.com/joey711/phyloseq/files/1795523/QIIME2_to_Phyloseq.pdf
# for converting qiime files into phyloseq format
metadata_ps <- read.csv(file="SeqAnalysis_mod/sample-metadata_phyloseq.csv", row.names = 1)

# Import taxonomy and count tables
count_table <- read_csv(file="asv_table.csv", skip = 1)
taxonomy <- read_csv(file="taxonomy_table.csv")
taxonomy <- taxonomy[-c(1),] 

# Remove characters and separate taxa names
taxonomy <-  taxonomy %>%
  mutate(taxonomy=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species"), sep=";")
colnames(taxonomy)[1] <- colnames(count_table)[1]

taxonomy_ps	=	as.matrix(taxonomy)
taxonomy_ps <- subset(taxonomy_ps, select = -c(Taxon, Confidence))
rownames(taxonomy_ps) <- subset(taxonomy_ps, select = c('#OTU ID'))
taxonomy_ps <- subset(taxonomy_ps, select = -c(1))


count_table_ps <- read.csv(file="SeqAnalysis_mod/table.csv", row.names = 1)
otu_table	=	as.matrix(count_table_ps)

phy_tree = read_tree("SeqAnalysis_mod/tree.nwk")

#	import	as	phyloseq	objects (phy_tree is already)
OTU	=	otu_table(count_table_ps,	taxa_are_rows	=	TRUE)
TAX	=	tax_table(taxonomy_ps)
META	=	sample_data(metadata_ps)

#	check	that	your	OTU	names	are	consistent	across	objects
taxa_names(TAX)
taxa_names(OTU)
taxa_names(phy_tree)

#	make	sure	files	have	the	same	sample	names	
sample_names(OTU)
sample_names(META)

#	merge	into	one	phyloseq	object
ps	=	phyloseq(OTU,	TAX,	META, phy_tree)
ps


# Now returning to using instructions from https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html
# Taxonomic Filtering
# Show available ranks in the dataset
rank_names(ps)

# Filter out those ambigious Kingdom annotations
ps <- subset_taxa(ps, !is.na(Kingdom) & !Kingdom %in% c("", "Unassigned"))

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Also filter out any with "NA" as phylum
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "Unassigned"))

# Check
unique(tax_table(ps)[, "Phylum"])




###### Eplore dataset a bit ##########

# Compute prevalence of each taxon, store as data.frame
# (prevalence = the number of samples in which a taxon appears at least once)
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})


# # NOTE- this next line takes a few minutes! uncomment only if need to rerun from scratch
# # Add taxonomy and total read counts to this data.frame
# prevdf = data.frame(Prevalence = prevdf, TotalAbundance = taxa_sums(ps), tax_table(ps))


# Are there phyla that are comprised of mostly low-prevalence features? 
# Compute the total and average prevalences of the features in each phylum.
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# Yes- in my case AncK6, Asgardaeota, others...
# In the tutorial they filter these out
# I am going to skip this for now because it might be interesting to see these differences in prevalence between sample types
# for some of these low abundance taxa


# explore the relationship of prevalence and total read count for each feature
p <- ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
ggsave("Figures/Prevalence_by_phylum.pdf",plot = p, dpi = 300, width = 12, height = 12, units = c("in"))
# This plots prevalence for every unique ASV (each dot) and divides the plot by Phylum
# This shows that there are some phyla (like Anck6) for which there are very few ASVs and furthermore,
# for some of them the ASVs only appear in a small percentage of samples
# There is a line in this plot as a cutoff of 5% for minimum prevalence
# Again, I am going to proceed without doing this but can come back to the tutorial 
# if I feel like this filtering is important.
# I suspect these low prevalence ASVs won't be important in the downstream analysis
# but in case they are, I can leave them in and check which phyla they are from later


# Abundance value transformation
# Skip ahead in tutorial
# This is for transforming abundance data to account for library size differences

# First, define plot_abundance(), which  uses phyloseq’s function to define a relative abundance graphic
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # This is an arbitrary subset, based on the Phylum Firmicutes, for plotting
  # (when doing this for Proteobacteria it takes a very long time)
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "Size_fraction",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform phyloseq object's otu counts to relative abundance. Save as new object.
ps_ra = transform_sample_counts(ps, function(x){x / sum(x)})

# Now we plot the abundance values before and after transformation.
plotBefore = plot_abundance(ps,"")
plotAfter = plot_abundance(ps_ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)

# Examine this figure - top panel shows original abundance without transformation. Below shows after transformation
# In the case of Firmicutes, it is not too different, suggesting that each sample reached saturation



# A more comprehensive way to do this is probably to check the rarefaction plots.
# Following this tutorial, https://micca.readthedocs.io/en/latest/phyloseq.html, 
# which recomends using the function from vegan

rarecurve(t(otu_table(ps)), step=50, cex=0.5)
# also takes a little bit of time to plot

# saved this figure quickly using GUI
# shows that all samples are saturated, even A1 which had lowest abundance of ASVs, seems to have just reached its plateau
# But still, keep an eye on this one to see if its ever an outlier or if composition looks similar to A2


# Check abundance  of all ASVs across samples to make sure there are no outliers in terms of low abundance samples
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
# distribution looks good




############## Ordinations ##############
# Try different ordinations to see which explains most variability

# First need to transform to reduce influence of zeroes using Hellinger transformation
# (= sqrt of rel abun)
ps_ra = transform_sample_counts(ps, function(x){(x)/sum(x)})
ps_hellinger <- transform_sample_counts(ps_ra, function(x){sqrt(x)})

# check
head(otu_table(ps), 25)
head(otu_table(ps_ra), 25)
head(otu_table(ps_hellinger), 25)



# PCoA (MDS)
out.pcoa.log <- ordinate(ps_hellinger, method = "MDS", distance = "bray")
evals <- out.pcoa.log$values$Eigenvalues
pcoa_plot = plot_ordination(ps_hellinger, out.pcoa.log, shape="Size_fraction") +
  geom_point(size = 3) +
  geom_text(mapping = aes(label = YYYYMMDD), size = 2)
            #, position = position_jitter(width = .05, height = .03))

ggsave(filename = "Figures/pcoa_plot.jpg", plot = pcoa_plot, dpi = 300, width = 4, height = 3, units = c("in"))

# Date labels show that replicates are very close to each other- good.
# Also the Pfree and Total are not far from each other

# Axis 1- 23.9%, axis2- 17.6%





# PcOA using weighted unifrac
out.wuf.log <- ordinate(ps_hellinger, method = "MDS", distance = "wunifrac")
evals <- out.wuf.log$values$Eigenvalues
wuf_plot = plot_ordination(ps_hellinger, out.wuf.log, shape="Size_fraction") +
  geom_point(size = 3) +
  geom_text(mapping = aes(label = YYYYMMDD), size = 2)
            #, position = position_jitter(width = .01, height = .01))

ggsave(filename = "Figures/wuf_plot.jpg", plot = wuf_plot, dpi = 300, width = 4, height = 3, units = c("in"))
# Axis 1- 33.7%, axis2- 20.5%

# double principal coordinates analysis (DPCoA)
# out.dpcoa.log <- ordinate(ps_hellinger, method = "DPCoA")
# evals <- out.dpcoa.log$eig
# plot_ordination(ps_hellinger, out.wuf.log) +
#   coord_fixed(sqrt(evals[2] / evals[1]))
# this took >1hr to run then gave an error. skip


# nMDS using weighted unifrac
out.nmds.log <- ordinate(ps_hellinger, method = "NMDS", distance = "wunifrac")
evals <- out.nmds.log$values$Eigenvalues
nmds_plot = plot_ordination(ps_hellinger, out.nmds.log, shape="Size_fraction") +
  geom_point(size = 3) +
  geom_text(mapping = aes(label = YYYYMMDD), size = 2)
#, position = position_jitter(width = .01, height = .01))

ggsave(filename = "Figures/nmds_plot.jpg", plot = nmds_plot, dpi = 300, width = 4, height = 3, units = c("in"))
# stress is 0.130




# I am going to continue with the weighted unifrac PCoA. Shows strongest axes (>60% in 1st two axes)




# Plot environmental variables on PCoA

# First Salinity- 2 ways of doing this

# See distribution frequency of salinity values
qplot(sample_data(ps_hellinger)$Salinity, geom = "histogram",binwidth=1) + xlab("Salinity")

# 1st option- Divide into two binds for plotting
sample_data(ps_hellinger)$salinity_binned <- cut(sample_data(ps_hellinger)$Salinity, breaks = c(0, 5, 10))
levels(sample_data(ps_hellinger)$salinity_binned) <- list(Fresh_0_5="(0,5]", Brackish_5_10="(5,10]")

plot_ordination(ps_hellinger, out.wuf.log, color = "salinity_binned") +
  labs(col = "Binned Salinity")


# 2nd option= plot as continuous scale bar. This is probably more fair
wuf_sal_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Salinity") +
  geom_point(size = 4) +
  labs(col = "Salinity") +
  scale_color_viridis(option="magma") 
ggsave(filename = "Figures/wuf_sal_plot.jpg", plot = wuf_sal_plot, dpi = 300, width = 4, height = 3, units = c("in"))



# Temperature
wuf_temp_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Temperature") +
  geom_point(size = 4) +
  labs(col = "Temp (C)") +
  scale_color_viridis(option="magma") 
ggsave(filename = "Figures/wuf_temp_plot.jpg", plot = wuf_temp_plot, dpi = 300, width = 4, height = 3, units = c("in"))

# axis 2 seems associated with temperature


# Secchi Depth
wuf_secchi_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "SecchiDepth") +
  geom_point(size = 4) +
  labs(col = "Secchi Dep (m)") +
  scale_color_viridis(option="magma") 
ggsave(filename = "Figures/wuf_secchi_plot.jpg", plot = wuf_secchi_plot, dpi = 300, width = 4, height = 3, units = c("in"))

# Note- not all points have a secchi measurement
# I don't see a pattern


# O2_mgL
wuf_O2_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "O2_mgL") +
  geom_point(size = 4) +
  labs(col = "O2 (mg/L)") +
  scale_color_viridis(option="magma") 
ggsave(filename = "Figures/wuf_O2_plot.jpg", plot = wuf_O2_plot, dpi = 300, width = 4, height = 3, units = c("in"))
# Axis 2 is sort of associated with O2, similar to temperature


# Volume filtered- see if this has an effect
wuf_volfilt_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Vol_Filtered") +
  geom_point(size = 4) +
  labs(col = "Vol.Filtered(ml)") +
  scale_color_viridis(option="magma") 
ggsave(filename = "Figures/wuf_volfilt_plot.jpg", plot = wuf_volfilt_plot, dpi = 300, width = 4, height = 3, units = c("in"))
# Volume filtered has no apparent effect so that's good!


# also color code by already existing categoriacal variables

# Location
wuf_location_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Location") +
  labs(col = "Location") +
  geom_point(size = 4)
ggsave(filename = "Figures/wuf_location_plot.jpg", plot = wuf_location_plot, dpi = 300, width = 4, height = 3, units = c("in"))


#Time of day
wuf_time_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Time") +
  labs(col = "Time of Day") +
  geom_point(size = 4)
ggsave(filename = "Figures/wuf_time_plot.jpg", plot = wuf_time_plot, dpi = 300, width = 4, height = 3, units = c("in"))


# Tide
wuf_tide_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Tide") +
  labs(col = "Tide") +
  geom_point(size = 4)
ggsave(filename = "Figures/wuf_tide_plot.jpg", plot = wuf_tide_plot, dpi = 300, width = 4, height = 3, units = c("in"))


#Rain
wuf_rain_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Rain_3d_prior") +
  labs(col = "Rain 3d Prior?") +
  geom_point(size = 4)
ggsave(filename = "Figures/wuf_rain_plot.jpg", plot = wuf_rain_plot, dpi = 300, width = 4, height = 3, units = c("in"))


#Size fraction
wuf_size_fraction_plot <- plot_ordination(ps_hellinger, out.wuf.log, color = "Size_fraction") +
  labs(col = "Size Fraction") +
  geom_point(size = 4)
ggsave(filename = "Figures/wuf_size_fraction_plot.jpg", plot = wuf_size_fraction_plot, dpi = 300, width = 4, height = 3, units = c("in"))
# This seems associated with Axis 1



# Differential Abundance
## There are many more useful things in this tutorial I may want to come back to
# But why I really wanted to phyloseq was to do Differential Abundance analyses

# Following tutorial from here https://joey711.github.io/phyloseq-extensions/DESeq2.html
# and https://micca.readthedocs.io/en/latest/phyloseq.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

# library("DESeq2")
# packageVersion("DESeq2") #1.26.0



# Convert to deseq
# For some reason I can only do this with the original count matrix and not the transformed, because
# it wants to only take integers.
# But this makes me very skeptical of results (comparing apples to oranges)

# What is better is to use rel abun but multiply the ps_ra values by 1000000 so everything is out of 1000000 instead of 1 
# That way, deseq will round it to integers (I think) and you'rr comparing normalize samples, 
# not samples with different sequencing effort
ps_ra_1000000 <- transform_sample_counts(ps_ra, function(x){x*1000000})
# check
head(otu_table(ps_ra), 25)
head(otu_table(ps_ra_1000000), 25)

#convert to deseq- using size fraction as the comparison variable
Size_fraction_dds = phyloseq_to_deseq2(ps_ra_1000000, ~Size_fraction)
Size_fraction_dds = DESeq(Size_fraction_dds, test="Wald", fitType="parametric")

# Investigate results table
res = results(Size_fraction_dds, cooksCutoff = FALSE)
# Filter out those taxa with pvals<0.01
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Now check out which ASVs were significantly different between sample types
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Look at genera only in Proteobacteria
sigtab_proteo <- subset(sigtab, Phylum %in% c("Proteobacteria"))
ggplot(sigtab_proteo, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Look at genera only in Firmicutes
sigtab_firmi <- subset(sigtab, Phylum %in% c("Firmicutes"))
ggplot(sigtab_firmi, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Look at genera only in Bacteroidetes
sigtab_bacteroid <- subset(sigtab, Phylum %in% c("Bacteroidetes"))
ggplot(sigtab_bacteroid, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Look at genera only in Actinobacteria
sigtab_actino <- subset(sigtab, Phylum %in% c("Actinobacteria"))
ggplot(sigtab_actino, aes(x=Genus, y=log2FoldChange, color=Order)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))





# Continue examining above 

# Also do analysis for rain variable (rain 3d prior? yes or no)

#convert to deseq- using size fraction as the comparison variable
# First filter out samples with rain = NA
Rain_dds = phyloseq_to_deseq2(ps_ra_1000000, ~Rain_3d_prior)
Rain_dds = DESeq(Rain_dds, test="Wald", fitType="parametric")
#### STOPPED HERE. FIX ABOVE. FILTER OUT NAs







# Continue here- assess the ASVs that are differentially distributed in Pfree vs Total
# Check the results from this figure
# Also reconsider using log2 fold change vs something else (that's really for expression data)
# Also what does positive cs negative mean here in terms of Pfree or Total
# Not seeing a lot of negatives in plot









### Bar Plots
# follow https://github.com/joey711/phyloseq/issues/901

# df <- psmelt(ps_ra)  # this line takes a while to run

# Define the number of colors you want
library(RColorBrewer)
nb.cols <- length(unique(df$Phylum))

# Build plot
p <- ggplot(data=df, aes(x=Sample, y=Abundance, fill=Phylum)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  theme(legend.position="bottom") + guides(fill=guide_legend(nrow=6)) +
  theme(legend.text = element_text(size = 8)) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)) +
  ylim(0,1)
p
ggsave(filename = "Figures/plylum_barplot.jpg", plot = p, dpi = 1200, width = 20, height = 10, units = c("in"))


### Bar Plots- alternative
# follow https://micca.readthedocs.io/en/latest/phyloseq.html




## Next
# Filter by thresholds to clean up plots- https://github.com/joey711/phyloseq/issues/901
# Clean up script
# Think about differential abund more
