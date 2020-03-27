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


## Load packages 

# library(stats)
# library(tidyverse)
# library(vegan)
# library(ggthemes)
# library(reshape2)
# library(devtools)



# rm(list=ls(all=TRUE))


# Import taxonomy and count tables
count_table <- read_csv(file="asv_table.csv", skip = 1)
taxonomy <- read_csv(file="taxonomy_table.csv")
taxonomy <- taxonomy[-c(1),] 

# Remove characters and separate taxa names
taxonomy <-  taxonomy %>%
  mutate(taxonomy=str_replace_all(string=Taxon, pattern="D_\\d*\\__", replacement="")) %>%
  mutate(taxonomy=str_replace_all(string=taxonomy, pattern=";$", replacement="")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus","species"), sep=";")
colnames(taxonomy)[1] <- colnames(count_table)[1]

# Join taxonomy to counts
colnames(taxonomy)[1] <- colnames(count_table)[1]
asv_count_table <- inner_join(taxonomy, count_table, by = colnames(taxonomy)[1])


# Calculate relative abundances based on total ASVs in each sample
samplesums <- colSums(asv_count_table[,c(11:63)],na.rm = TRUE)
samplesums
for(i in 1:ncol(asv_count_table[,c(11:63)])){
  asv_count_table[,i+10] <- asv_count_table[,i+10]/samplesums[i]
}



# Calculate averages for replicate samples. Put into new table
asv_relabun <- asv_count_table
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(A=mean(c(A1, A2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(B=mean(c(B1, B2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(C=mean(c(C1, C2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(D=mean(c(D1, D2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(E=mean(c(E1, E2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(F=mean(c(F1, F2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(G=mean(c(G1, G2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(H=mean(c(H1, H2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(I=mean(c(I1, I2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(J=mean(c(J1, J2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(K=mean(c(K1, K2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(L=mean(c(L1, L2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(M=mean(c(M1, M2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(N=mean(c(N1, N2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(O=mean(c(O1, O2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(P=mean(c(P1, P2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(Q=mean(c(Q1, Q2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(R=mean(c(R1, R2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(S=mean(c(S1, S2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(T=mean(c(T1, T2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(U=mean(c(U1, U2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(V=mean(c(V1, V2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(W=mean(c(W1, W2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(X=mean(c(X1, X2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(Y=mean(c(Y1, Y2), na.rm=T)) 
asv_relabun <- asv_relabun %>% rowwise() %>% mutate(Z=mean(c(Z1, Z2), na.rm=T)) 

asv_relabun <- select(asv_relabun, -A1, -A2, -B1, -B2, -C1, -C2, -D1, -D2, -E1, -E2, -F1, -F2, -G1, -G2, -H1, -H2, -I1, -I2, -J1, -J2, -K1, -K2, -L1, -L2, -M1, -M2, -N1, -N2, -O1, -O2, -P1, -P2, -Q1, -Q2, -R1, -R2, -S1, -S2, -T1, -T2, -U1, -U2, -V1, -V2, -W1, -W2, -X1, -X2, -Y1, -Y2, -Z1, -Z2)


# Make preliminary bubble plots

# Phyla
phyla <- unique(asv_count_table$phylum)
phyla
# There are 53 unique phyla- that's a lot! But the annotation is not always great so could be some error in this. 
# Check carefully 
# remove "NA" phylum so there are now only 52 phyla
phyla <- phyla[-c(11)]
# Make an empty matrix
phyla_table <-  matrix(, nrow = length(phyla), ncol = 26) 
rownames(phyla_table) <- phyla
colnames(phyla_table) <- colnames(asv_relabun[,12:37])
# Fill empty matrix with sums of that respective phylum and sample 
for(i in 1:length(phyla)){
  x <- which(asv_relabun$phylum == phyla[i]) # x is index (rows) of all instances of that taxon
  for(y in 12:37){
    phyla_table[i,y-11] <- sum(asv_relabun[x,y], na.rm = TRUE)
  }
}
# Convert to tibble (keeping row names with workaround)
phyla_table <- as_tibble(phyla_table, rownames = "phyla")
# Reshape the matrix to make it suitable for ggplot
phyla_table <- melt(phyla_table, id.vars = "phyla", variable.name="Samples", value.name = "RelAbun")
#Turn your phyla column into a character vector
phyla_table$phyla <- as.character(phyla_table$phyla)
#Then turn it back into a factor with the levels in the correct order
phyla_table$phyla <- factor(phyla_table$phyla, levels=sort(unique(phyla_table$phyla)))
# Plot with factors (phyla) in 'reverse' order so they y=axis appears alphabetically top to btm
phyla_plot <- ggplot(phyla_table, aes(x = Samples, y = fct_rev(phyla))) +  
  geom_point(aes(size = RelAbun)) + 
  scale_size(range = range(c(-1, 2)),limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Phyla")
phyla_plot
ggsave("Figures/phyla_plot.jpg", plot = phyla_plot)




# Class
class <- unique(asv_count_table$class)
class
# There are 142 unique "classes"
# remove "NA" 
class <- class[!is.na(class)]
# Make an empty matrix
class_table <-  matrix(, nrow = length(class), ncol = 26) 
rownames(class_table) <- class
colnames(class_table) <- colnames(asv_relabun[,12:37])
# Fill empty matrix with sums of that respective phylum and sample 
for(i in 1:length(class)){
  x <- which(asv_relabun$class == class[i]) # x is index (rows) of all instances of that taxon
  for(y in 12:37){
    class_table[i,y-11] <- sum(asv_relabun[x,y], na.rm = TRUE)
  }
}
# Convert to tibble (keeping row names with workaround)
class_table <- as_tibble(class_table, rownames = "class")
# Reshape the matrix to make it suitable for ggplot
class_table <- melt(class_table, id.vars = "class", variable.name="Samples", value.name = "RelAbun")
#Turn your class column into a character vector
class_table$class <- as.character(class_table$class)
#Then turn it back into a factor with the levels in the correct order
class_table$class <- factor(class_table$class, levels=sort(unique(class_table$class)))
# Plot with factors (class) in 'reverse' order so they y=axis appears alphabetically top to btm
class_plot <- ggplot(class_table, aes(x = Samples, y = fct_rev(class))) +  
  geom_point(aes(size = RelAbun)) + 
  scale_size(range = range(c(-1, 2)),limits = c(0,1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Classes")
class_plot
ggsave("Figures/class_plot.jpg", plot =class_plot, width = 13, height = 13, units = "in")



## TO DO: remove low abundance phyla and classes

### NEXT- PLOT BARPLOTS IN COLOR (?) FOR ELISE


# Load metadata
metadata <- read_csv(file="sample-metadata-R.csv")

# Subset and organize asv_relabun so samples are rows and ASV IDs are columns
# First subset
asv_table <- asv_relabun %>%
  select('#OTU ID', A, B, C, D, E, F, G, H, I ,J, K, L, M ,N, O, P, Q, R, S, T, U, V, W, X, Y, Z)
# convert to dataframe (necessary for transformation)
asv_table <- as.data.frame(asv_table)
# give ASVs names as row names 
rownames(asv_table) <- asv_table[,1]
asv_table <- asv_table[-c(1)]
# and transform
asv_table <- t(asv_table)


# Calculate nMDS
NMDS_asv_table=metaMDS(asv_table)
# Stress is ~0.15. This is acceptable
plot(NMDS_asv_table, display = "sites", xlim=c(-2.1, 1), ylim=c(-2, 2),cex.lab = 1, cex.axis = 1) # displays by samples

# Do env fit of metadata to nmds
attach(metadata)
NMDS_asv_table_envfit <- envfit(NMDS_asv_table~Year+Month+Day+YYYYMMDD+Location+Salinity+Temperature+Vol_Filtered++Size_fraction, permutations = 9999)
# Left out Time+SecchiDepth+Tide+O2_%+O2_mgL+Rain_3d_prior because 1) Can't handle number format of time and 2) Can't handle NAs (bring up with Elise- get these data?)

# # Save stats results
# capture.output(NMDS_asv_table, file = "StatResults/NMDS_asv_table.txt")
# capture.output(NMDS_asv_table_envfit, file = "StatResults/NMDS_asv_table_envfit.txt")
# # Salinity, Temp, Month are significant

# Plot "sites" color-coded by Size fraction and labelled by date
bitmap("Figures/nmds_plot.jpg", height = 4, width = 4, units = 'in', type="tifflzw", res=1200)
with(metadata, levels(as.factor(Size_fraction)))
scl <- 2
colvec <- c("black","lightgray") # Pfree is black, total is grey
plot(NMDS_asv_table, display = "sites", xlim=c(-2.1, 1), ylim=c(-2, 2),cex.lab = 1, cex.axis = 1) # displays by samples
box(lwd=2)
text(NMDS_asv_table, labels=YYYYMMDD, pos = 2, cex=0.9)
with(metadata, points(NMDS_asv_table, display = "sites", col = colvec[as.factor(Size_fraction)], pch = 21, bg = colvec[as.factor(Size_fraction)]), cex=1.5)
legend(.6, -1.3, legend=levels(as.factor(Size_fraction)), col=colvec, pch = 19, cex = 0.7)
dev.off()
par(mfrow = c(1,1))


### NEXT- PLOT PCOA FOR ELISE






### Plot abundance of taxa against each other for paired Total/Pfree samples to see which taxa are not represented the same in each sample type

# 20180226 - sample B (Pfree) and S (Total)- These are very different on the nMDS
# Plot with factors (class) in 'reverse' order so they y=axis appears alphabetically top to btm
ggplot(asv_relabun, aes(x = B, y = S)) +  
  geom_point() 

# 20170613- sample H (total) vs L (Pfree)- Look very similar on nMDS
ggplot(asv_relabun, aes(x = H, y = L)) +  
  geom_point() 
# These look as expected. Now try to figure out which taxa are the outliers in each case

# Calculate a ratio, for each taxon, of Pfree:Total and set a threshold (eg. if >2 or <.5)???
temp <- asv_relabun$S/asv_relabun$B
temp[which(temp == Inf)] <- NA
# temp[which(temp == NaN)] <- NA

asv_relabun$Taxon[c(which(temp > 2))]
##CONTINUE HERE



# ALTERNATIVE-
# try microbiomeseq package 
# first need to convert data to phyloseq format
# to use phyloseq (https://github.com/joey711/phyloseq)
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install_github("umerijaz/microbiomeSeq")  
# install.packages("microbiomeSeq")
# library(microbiomeSeq) 

# NOTE. I can't get this to install. Developer says this is still beta and they don't recommend using it until they give it OK on this page:
# http://userweb.eng.gla.ac.uk/umer.ijaz/projects/microbiomeSeq_Tutorial.html#content



# Another option is from phyloseq but also requires changing file formats
# https://joey711.github.io/phyloseq-extensions/DESeq2.html
# OR another option is permubiome
# https://cran.r-project.org/web/packages/permubiome/permubiome.pdf
















# Plot abundances of target organisms genera: Enterococcus, Salmonella, Shigella, Vibrio, Actinomycetes

# Entero
x <- which(asv_relabun$genus == "Enterococcus") 
x
# Row 6159 is only row with enterococcus. Checked the ID by blasting this sequence (ID: 88bad5fd626727f0d8fef57f9dbb6e0f) and it is in the Enterococcus genus
# Saved blast hits in /BlastResults as enterococcus.csv
entero <- asv_relabun[x,9:37]
entero <- entero[,-c(2,3)]
entero <- melt(entero, id.vars = "genus", variable.name="Samples", value.name = "RelAbun")
colnames(metadata)[1] <- colnames(entero)[2]
entero <- left_join(entero, metadata, by = "Samples")
# Plot by sample ID
ggplot(entero, aes(x = Samples, y = RelAbun)) +  
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# Plot by month- PLAY AROUND HERE WITH MONTH/ YEAR ETC
ggplot(entero, aes(x = Month, y = RelAbun)) +  
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Salmonella
x <- which(asv_relabun$genus == "Salmonella") 
x
# Doesn't exist


# Shigella
x <- which(asv_relabun == "Shigella") 
x
# Doesn't exist


# Vibrio
x <- which(asv_relabun$genus == "Vibrio") 
x
Vibrio <- asv_relabun[x,]
# 10 rows have Vibrio (and two of these are Vibrio Cholerae)
# Make a Vibrio plot and a Vibrio cholerae plot
vibrio <- asv_relabun[x,9:37]
vibrio <- vibrio[,-c(2,3)]

# Sum the abundance of all Vibrio for each sample:
# First make empty vector
vibrio_table <-  matrix(, nrow = 1, ncol = 26) 
rownames(vibrio_table) <- "Vibrio"
colnames(vibrio_table) <- colnames(asv_relabun[,12:37])
# Fill empty matrix with sums of that respective phylum and sample 
for(y in 2:27){
  vibrio_table[1,y-1] <- sum(vibrio[,y], na.rm = TRUE)
  }
# Plot
vibrio_table <- melt(vibrio_table, id.vars = "genus", variable.name="Samples", value.name = "RelAbun")
colnames(vibrio_table)[1] <- "genus"
colnames(vibrio_table)[2] <- "Samples"
ggplot(vibrio_table, aes(x = Samples, y = RelAbun)) +  
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Now sum the abundance of just Vibrio Cholerae for each sample:
x <- which(asv_relabun$species == "Vibrio cholerae") 
x
asv_relabun[x,]
# row 353 (id 9fe7b6afa4d0f3a3347aba994cffa125) and 1556 (id 6564a6ce8fa69294e0a9bfb5c97b717f)
# saved results in /BlastResults as Virbrio1.csv and Vibrio2.csv respectively
# note that in blast, V. cholerae and V. albensis (as well as some other Vibrio) are 100% identical in this V4 region
vibrio_chol <- asv_relabun[x,9:37]
vibrio_chol <- vibrio_chol[,-c(1,3)]
#  Make empty vector
vibrio_chol_table <-  matrix(, nrow = 1, ncol = 26) 
rownames(vibrio_chol_table) <- "Vibrio_cholerae"
colnames(vibrio_chol_table) <- colnames(asv_relabun[,12:37])
# Fill empty matrix with sums of that respective phylum and sample 
for(y in 2:27){
  vibrio_chol_table[1,y-1] <- sum(vibrio_chol[,y], na.rm = TRUE)
}
# Plot
vibrio_chol_table <- melt(vibrio_chol_table, id.vars = "species", variable.name="Samples", value.name = "RelAbun")
colnames(vibrio_chol_table)[1] <- "species"
colnames(vibrio_chol_table)[2] <- "Samples"
ggplot(vibrio_chol_table, aes(x = Samples, y = RelAbun)) +  
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))











