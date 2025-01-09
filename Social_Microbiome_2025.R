#### Analysis of sociality and gut microbiome in California ground squirrels #####
### by Erin Person, January 2025 ###
### Adapted from Callahan et al. 2017 on Bioconductor
### https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html

#### SET-UP ####

# clear workspace
rm(list = ls())
# set seed
set.seed(42)

# load packages
library(dada2)
library(phyloseq)
library(phyloseqCompanion)
library(Biostrings)
library(tidyverse)
library(lme4)
library(DECIPHER)
library(phangorn)
library(decontam)
library(DESeq2)
library(lmerTest)
library(asnipe)
library(igraph)
library(ANTs)

# set working directory and path
setwd("your path here")
path <- "your path here"

#### IMPORT DATA ####

# Reads in file names following the format described in the string
# This line creates a list of the forward reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#### PROCESSING IN DADA2 ####

# plot read quality profiles for forward reads for the first two samples
plotQualityProfile(fnFs[1:2])

# plot read quality profiles for reverse reads
plotQualityProfile(fnRs[1:2])

# Create new filtered subdirectory and place filtered files there
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# truncating and filtering reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,150),
                     maxN=0, maxEE=c(3,5), truncQ=2, trimLeft = c(10,10), rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

# dereplication (reduces computation time)
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# name the derep-class objects by sample names
names(derepFs) <- sample.names
names(derepFs) <- sample.names

# learning error rates (forward and reverse)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# plot errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# applies core sample inference algorithm to trimmed and filtered data
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = FALSE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = FALSE)

# merges forward and reverse reads into contigs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# construct amplicon sequence variant (ASV) table
seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#track reads through pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# assign taxonomy with Silva db
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/perse/Documents/Berkeley/Projects/Ground squirrel sociality and microbiome/Microbiome Chapter/P-00K9 MI Project/16Sv4_Sequences/silva_nr99_v138_train_set.fa.gz", 
                       multithread=TRUE, verbose = TRUE)

# inspect taxonomic assignments
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# creating a taxonomy tree in DECIPHER for use in phangorn tree
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # this propagates to the tip labels of the tree
alignment <- DECIPHER::AlignSeqs(Biostrings::DNAStringSet(seqs), anchor=NA,verbose=TRUE)

# use phangorn to construct a tree
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

#### PHYLOSEQ ####

#creating the phyloseq objects
OTU = otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX = tax_table(taxa)
PHY = phy_tree(fitGTR$tree)

#reading in metadata
#including rownames to match OTU table
samdf <- data.frame(read_csv('metadata_for_pub.csv'))
rownames(samdf) <- sample_names(OTU)

#finishing phyloseq objects
SAM = sample_data(samdf)
ps = phyloseq(OTU, TAX, SAM, PHY)

#creates short names for ASVs instead of long DNA strings
#DNA sequences can be found using "refseq" function (refseq(ps))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

### same as above but includes unID phyla
# Compute prevalence of each feature, store as data.frame
prevdf_unid = apply(X = otu_table(ps),
                    MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_unid = data.frame(Prevalence = prevdf_unid,
                         TotalAbundance = taxa_sums(ps),
                         tax_table(ps))

phylum_abun_filtered <- prevdf_unid %>%
  group_by(Phylum) %>%
  summarise(Abundance = sum(TotalAbundance))

#remove un-characterized phyla
ps1 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# decontam process
df <- as.data.frame(sample_data(ps1)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps1)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sample.group)) + geom_point()

# find contaminants at standard threshold (0.05)
sample_data(ps1)$is.neg <- sample_data(ps1)$sample.group == "control"
contamdf.prev <- isContaminant(ps1, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)

#  identify and remove contaminants
contaminants <- filter(contamdf.prev, contaminant == TRUE)
ps2 <- prune_taxa(!contamdf.prev$contaminant, ps1)

# identify mitochondrial sequences
mitochondrial.seqs <- subset_taxa(ps2, Family == "Mitochondria")
mitochondria.taxa <- colnames(otu_table(mitochondrial.seqs))

# identify non-bacterial sequences
nonbacterial.seqs <- subset_taxa(ps2, Kingdom == "Archaea")
nonbacterial.taxa <- colnames(otu_table(nonbacterial.seqs))

# identify chloroplast sequences
chloroplast.seqs <- subset_taxa(ps2, Order == "Chloroplast")
chloroplast.taxa <- colnames(otu_table(chloroplast.seqs))

#Function to remove bad taxa with thanks to Sondra Turjeman for code from
#"Comparing invasive and noninvasive faecal sampling in wildlife microbiome studies: A case study on wild common cranes"
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
#Apply function
#removed 6 mitochondrial ASVs, 5 chloroplast and 4 archaea/non-bacterial ASVs
ps2_filtered = pop_taxa(ps2,mitochondria.taxa)
ps2_filtered = pop_taxa(ps2_filtered,nonbacterial.taxa)
ps2_filtered = pop_taxa(ps2_filtered,chloroplast.taxa)

# remove negative controls
ps3 <- prune_samples(sample_data(ps2_filtered)$sample.group  == "sample", ps2_filtered)

# Compute prevalence of each feature, store as data.frame
prevdf_filtered = apply(X = otu_table(ps3),
                        MARGIN = ifelse(taxa_are_rows(ps3), yes = 1, no = 2),
                        FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_filtered = data.frame(Prevalence = prevdf_filtered,
                             TotalAbundance = taxa_sums(ps3),
                             tax_table(ps3))

phylum_abun_filtered <- prevdf_filtered %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(Abundance = sum(TotalAbundance))

table(tax_table(ps3)[, "Phylum"], exclude = NULL)

### rarefaction ###

# rarefy to 2000 reads
ps_2000 <- rarefy_even_depth(ps3, sample.size = 2000, rngseed = 42)

# Compute prevalence of each feature, store as data.frame
prevdf_2000 = apply(X = otu_table(ps_2000),
                    MARGIN = ifelse(taxa_are_rows(ps_2000), yes = 1, no = 2),
                    FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf_2000 = data.frame(Prevalence = prevdf_2000,
                         TotalAbundance = taxa_sums(ps_2000),
                         tax_table(ps_2000))

phylum_abun_rarefaction <- prevdf_2000 %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(Abundance = sum(TotalAbundance))

#### META DATA ####

# extract rarefied sample data from ps object
samdf_rare_all <- sample.data.frame(sample_data(ps_2000))
row.names(samdf_rare_all) <- samdf_rare_all$sample.id
samdf_rare_all$sample.date <- as.Date(samdf_rare_all$sample.date, format = "%e-%b-%y")
samdf_rare_all$year <- year(ymd(samdf_rare_all$sample.date))
samdf_rare_all$year <- as.factor(samdf_rare_all$year)
samdf_rare_all$month <- lubridate::month(ymd(samdf_rare_all$sample.date), label = TRUE)

# remove secondary site samples
samdf_rare <- samdf_rare_all %>%
  filter(site == "primary") 

# filtering to single sample per individual per year, closest to June 15
dupes <- samdf_rare %>%
  dplyr::group_by(year, uid) %>%
  dplyr::filter(n() > 1) %>% 
  dplyr::filter(!grepl("Jun",month)) %>%
  ungroup()

# filter to samples in June, now one sample/individual/year
samdf_rare <- subset(samdf_rare, !(sample.id %in% dupes$sample.id))

# creating shannon diversity column for rarefied samples only
diversity <- estimate_richness(ps_2000, split = TRUE, measures = c("Observed", "Shannon"))
diversity$sample.id <- row.names(diversity)
samdf_rare <- merge(samdf_rare, diversity, by = "sample.id")

# subsetting ps to df above
ps_2000A <- prune_samples(samdf_rare$sample.id, ps_2000)

# Calculating dissimilarity values
w_dist <- phyloseq::distance(ps_2000A, method = "unifrac", weighted=TRUE)
u_dist <- phyloseq::distance(ps_2000A, method = "unifrac", weighted=FALSE)

# ordinating dissimilarity values
w_ord <- ordinate(ps_2000A, method="PCoA", distance = w_dist)
u_ord <- ordinate(ps_2000A, method="PCoA", distance = u_dist)

# create axis columns
samdf_rare$Axis.1.w = w_ord$vectors[, 1]
samdf_rare$Axis.2.w = w_ord$vectors[, 2]
samdf_rare$Axis.1.u = u_ord$vectors[, 1]
samdf_rare$Axis.2.u = u_ord$vectors[, 2]

# create year subsets
samdf18 <- samdf_rare %>%
  dplyr::filter(year == "2018")
ps18 <- prune_samples(samdf18$sample.id, ps_2000)
u_dist18 <- phyloseq::distance(ps18, method = "unifrac", weighted=FALSE)
u18 <- ordinate(ps18, method="PCoA", distance = u_dist18)

w_dist18 <- phyloseq::distance(ps18, method = "unifrac", weighted=TRUE)
w18 <- ordinate(ps18, method="PCoA", distance = w_dist18)

samdf19 <- samdf_rare %>%
  dplyr::filter(year == "2019")
ps19 <- prune_samples(samdf19$sample.id, ps_2000)
u_dist19 <- phyloseq::distance(ps19, method = "unifrac", weighted=FALSE)
u19 <- ordinate(ps19, method="PCoA", distance = u_dist19)

w_dist19 <- phyloseq::distance(ps19, method = "unifrac", weighted=TRUE)
w19 <- ordinate(ps19, method="PCoA", distance = w_dist19)

samdf20 <- samdf_rare %>%
  dplyr::filter(year == "2020")
ps20 <- prune_samples(samdf20$sample.id, ps_2000)
u_dist20 <- phyloseq::distance(ps20, method = "unifrac", weighted=FALSE)
u20 <- ordinate(ps20, method="PCoA", distance = u_dist20)

w_dist20 <- phyloseq::distance(ps20, method = "unifrac", weighted=TRUE)
w20 <- ordinate(ps20, method="PCoA", distance = w_dist20)

samdf21 <- samdf_rare %>%
  dplyr::filter(year == "2021")
ps21 <- prune_samples(samdf21$sample.id, ps_2000)
u_dist21 <- phyloseq::distance(ps21, method = "unifrac", weighted=FALSE)
u21 <- ordinate(ps21, method="PCoA", distance = u_dist21)

w_dist21 <- phyloseq::distance(ps21, method = "unifrac", weighted=TRUE)
w21 <- ordinate(ps21, method="PCoA", distance = w_dist21)

samdf22 <- samdf_rare %>%
  dplyr::filter(year == "2022")
ps22 <- prune_samples(samdf22$sample.id, ps_2000)
u_dist22 <- phyloseq::distance(ps22, method = "unifrac", weighted=FALSE)
u22 <- ordinate(ps22, method="PCoA", distance = u_dist22)

w_dist22 <- phyloseq::distance(ps22, method = "unifrac", weighted=TRUE)
w22 <- ordinate(ps22, method="PCoA", distance = w_dist22)

# site subsets 
# create subsets of site data
  
# secondary site
samdf_secondary <- samdf_rare_all %>%
  dplyr::filter(site == "secondary") %>%
    dplyr::filter(uid != "R1173" & uid != "R1157") # remove two individuals that have been caught in primary site
    
# primary site  
samdf_primary <- samdf_rare_all %>%
  dplyr::filter(site == "primary")

# randomly subset primary site to N of secondary site
samdf_primary_subset <- samdf_primary[sample(nrow(samdf_primary), size=21), ]
samdf_site <- rbind(samdf_secondary, samdf_primary)
samdf_site_subset <- rbind(samdf_secondary, samdf_primary_subset)

# prune ps to site df
ps_2000B <- prune_samples(samdf_site$sample.id, ps_2000)
ps_2000C <- prune_samples(samdf_site_subset$sample.id, ps_2000)

# unifrac for all
# Calculate weighted UniFrac distances
w_dist_site <- phyloseq::distance(ps_2000B, method = "unifrac", weighted=TRUE)
w_site <- ordinate(ps_2000B, method="PCoA", distance = w_dist_site)
# create weighted axis columns
samdf_site$Axis.1.w = w_site$vectors[, 1]
samdf_site$Axis.2.w = w_site$vectors[, 2]

# Calculate unweighted UniFrac distances
u_dist_site <- phyloseq::distance(ps_2000B, method = "unifrac", weighted=FALSE)
u_site <- ordinate(ps_2000B, method="PCoA", distance = u_dist_site)
# create unweighted axis columns
samdf_site$Axis.1.u = u_site$vectors[, 1]
samdf_site$Axis.2.u = u_site$vectors[, 2]

# subsetted Unifrac
# Calculate weighted UniFrac distances
w_dist_site_subset <- phyloseq::distance(ps_2000C, method = "unifrac", weighted=TRUE)
w_site_subset <- ordinate(ps_2000C, method="PCoA", distance = w_dist_site_subset)
# create weighted axis columns
samdf_site_subset$Axis.1.w = w_site_subset$vectors[, 1]
samdf_site_subset$Axis.2.w = w_site_subset$vectors[, 2]

# Calculate unweighted UniFrac distances
u_dist_site_subset <- phyloseq::distance(ps_2000C, method = "unifrac", weighted=FALSE)
u_site_subset <- ordinate(ps_2000C, method="PCoA", distance = u_dist_site_subset)
# create unweighted axis columns
samdf_site_subset$Axis.1.u = u_site_subset$vectors[, 1]
samdf_site_subset$Axis.2.u = u_site_subset$vectors[, 2]

# alpha diversity
diversity_site <- estimate_richness(ps_2000B, split = TRUE, measures = c("Observed", "Shannon"))
diversity_site$sample.id <- row.names(diversity_site)
samdf_site <- merge(samdf_site, diversity_site, by = "sample.id")

#### STATS ####

# function to test for random effects
# likelihood ratio test between lm and lmer tests
lrt=function (m, m1) {
  L0=logLik(m)
  L1=logLik(m1)
  L01=as.vector(- 2 * (L0 - L1))
  df=attr(L1, "df") - attr(L0, "df")
  list(L01 = L01, df = df,
       "p-value" = pchisq(L01, df, lower.tail = FALSE))}

# linear models

# NOTE: some individuals are missing social data, these are automatically excluded
# from linear models with social stats
# scale variables
samdf_scale <- transform(samdf_rare,
                        degree.s = scale(degree),
                        betweenness.s = scale(betweenness),
                        eigen.s = scale(eigen),
                        stan.strength.s = scale(stan.strength),
                        strength.s = scale(strength),
                        range.area.s = scale(range.area))

# betweenness
m1.b <- lmer(Shannon ~ betweenness.s + range.area.s + age.class + sex + year + (1|uid), data = samdf_scale, REML = FALSE)
summary(m1.b)

# uid effect
m1.b_no_uid <- lm(Shannon ~ betweenness.s + range.area.s + age.class + sex + year, data = samdf_scale)
lrt(m1.b_no_uid, m1.b)

# degree
m1.d <- lmer(Shannon ~ degree.s + range.area.s + age.class + sex + year + (1|uid), data = samdf_scale, REML = FALSE)
summary(m1.d)

# uid effect
m1.d_no_uid <- lm(Shannon ~ degree.s + range.area.s + age.class + sex + year, data = samdf_scale)
lrt(m1.d_no_uid, m1.d)

# eigen
m1.e <- lmer(Shannon ~ eigen.s + range.area.s + age.class + sex + year + (1|uid), data = samdf_scale, REML = FALSE)
summary(m1.e)

# uid effect
m1.e_no_uid <- lm(Shannon ~ eigen.s + range.area.s + age.class + sex + year, data = samdf_scale)
lrt(m1.e_no_uid, m1.e)

# strength
m1.s <- lmer(Shannon ~ stan.strength.s + range.area.s + age.class + sex + year + (1|uid), data = samdf_scale, REML = FALSE)
summary(m1.s)

# uid effect
m1.s_no_uid <- lm(Shannon ~ stan.strength.s + range.area.s + age.class + sex + year, data = samdf_scale)
lrt(m1.s_no_uid, m1.s)

# site
m2 <- lm(Shannon ~ site, data = samdf_site)
summary(m2)

# permanovas

# Adonis test
p1 <- adonis2(u_dist ~ uid + year + social.group + age.class, data = samdf_rare, na.action = na.exclude, by = "terms", permutations = 1000)
p1

p2 <- adonis2(w_dist ~ uid + year + social.group + age.class, data = samdf_rare, na.action = na.exclude, by = "terms", permutations = 1000)
p2

# permanova by site
p3 <- adonis2(w_dist_site ~ site, by = "terms", data = samdf_site)
p3

p4 <- adonis2(u_dist_site ~ site, by = "terms", data = samdf_site)
p4

# subset
p5 <- adonis2(w_dist_site_subset ~ site, by = "terms", data = samdf_site_subset)
p5

p6 <- adonis2(u_dist_site_subset ~ site, by = "terms", data = samdf_site_subset)
p6

## comparing ai, hro and dissimilarity matrices

AI_phys <- read.csv("Phys_AI_2018-2022.csv")
AI_phys$year <- as.factor(AI_phys$year)

hro <- read.csv("percent_overlap_by_year.csv")
hro$year <- as.factor(hro$year)

# merging with samdf to match up sample ids in AI
AI2 <- left_join(AI_phys, samdf_rare, by = c("uid_ini" = "uid", "year"))
AI2 <- AI2 %>%
  dplyr::rename(sample.id.ini = sample.id) %>%
  dplyr::filter(!is.na(sample.id.ini))
AI3 <- left_join(AI2, samdf_rare, by = c("uid_rec" = "uid", "year"))
AI3 <- AI3 %>%
  dplyr::rename(sample.id.rec = sample.id) %>%
  dplyr::filter(!is.na(sample.id.rec)) %>%
  dplyr::select(uid_ini, sample.id.ini, uid_rec, sample.id.rec, AI, year)

AI18 <- AI3 %>%
  filter(year == "2018")
AIm18 <- df.to.mat(AI18, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "AI", sym = TRUE)
um18 <- as.matrix(u_dist18)
um18 <- um18[row.names(um18) %in% row.names(AIm18), colnames(um18) %in% colnames(AIm18)]
wm18 <- as.matrix(w_dist18)
wm18 <- wm18[row.names(wm18) %in% row.names(AIm18), colnames(wm18) %in% colnames(AIm18)]

AI19 <- AI3 %>%
  filter(year == "2019")
AIm19 <- df.to.mat(AI19, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "AI", sym = TRUE)
um19 <- as.matrix(u_dist19)
um19 <- um19[row.names(um19) %in% row.names(AIm19), colnames(um19) %in% colnames(AIm19)]
wm19 <- as.matrix(w_dist19)
wm19 <- wm19[row.names(wm19) %in% row.names(AIm19), colnames(wm19) %in% colnames(AIm19)]

AI20 <- AI3 %>%
  filter(year == "2020")
AIm20 <- df.to.mat(AI20, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "AI", sym = TRUE)
um20 <- as.matrix(u_dist20)
um20 <- um20[row.names(um20) %in% row.names(AIm20), colnames(um20) %in% colnames(AIm20)]
wm20 <- as.matrix(w_dist20)
wm20 <- wm20[row.names(wm20) %in% row.names(AIm20), colnames(wm20) %in% colnames(AIm20)]

AI21 <- AI3 %>%
  filter(year == "2021")
AIm21 <- df.to.mat(AI21, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "AI", sym = TRUE)
um21 <- as.matrix(u_dist21)
um21 <- um21[row.names(um21) %in% row.names(AIm21), colnames(um21) %in% colnames(AIm21)]
wm21 <- as.matrix(w_dist21)
wm21 <- wm21[row.names(wm21) %in% row.names(AIm21), colnames(wm21) %in% colnames(AIm21)]

AI22 <- AI3 %>%
  filter(year == "2022")
AIm22 <- df.to.mat(AI22, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "AI", sym = TRUE)
um22 <- as.matrix(u_dist22)
um22 <- um22[row.names(um22) %in% row.names(AIm22), colnames(um22) %in% colnames(AIm22)]
wm22 <- as.matrix(w_dist22)
wm22 <- wm22[row.names(wm22) %in% row.names(AIm22), colnames(wm22) %in% colnames(AIm22)]

# same as above but for home range overlap
hro2 <- left_join(hro, samdf_rare, by = c("Var1" = "uid", "year"))
hro2 <- hro2 %>%
  dplyr::rename(sample.id.ini = sample.id) %>%
  filter(!is.na(sample.id.ini))
hro3 <- left_join(hro2, samdf_rare, by = c("Var2" = "uid", "year"))
hro3 <- hro3 %>%
  dplyr::rename(sample.id.rec = sample.id) %>%
  filter(!is.na(sample.id.rec)) %>%
  dplyr::select(Var1, sample.id.ini, Var2, sample.id.rec, value, year) %>%
  dplyr::rename(uid.ini = Var1, uid.rec = Var2, overlap = value)

hro18 <- hro3 %>%
  filter(year == "2018")
hrom18 <- df.to.mat(hro18, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "overlap", sym = FALSE)
hrom18 <- hrom18[row.names(hrom18) %in% row.names(AIm18), colnames(hrom18) %in% colnames(AIm18)]

hro19 <- hro3 %>%
  filter(year == "2019")
hrom19 <- df.to.mat(hro19, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "overlap", sym = FALSE)
hrom19 <- hrom19[row.names(hrom19) %in% row.names(AIm19), colnames(hrom19) %in% colnames(AIm19)]

hro20 <- hro3 %>%
  filter(year == "2020")
hrom20 <- df.to.mat(hro20, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "overlap", sym = FALSE)
hrom20 <- hrom20[row.names(hrom20) %in% row.names(AIm20), colnames(hrom20) %in% colnames(AIm20)]

hro21 <- hro3 %>%
  filter(year == "2021")
hrom21 <- df.to.mat(hro21, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "overlap", sym = FALSE)
hrom21 <- hrom21[row.names(hrom21) %in% row.names(AIm21), colnames(hrom21) %in% colnames(AIm21)]

hro22 <- hro3 %>%
  filter(year == "2022")
hrom22 <- df.to.mat(hro22, actor = "sample.id.ini", receiver = "sample.id.rec", weighted = "overlap", sym = FALSE)
hrom22 <- hrom22[row.names(hrom22) %in% row.names(AIm22), colnames(hrom22) %in% colnames(AIm22)]

# prune association index and unifrac matrices to same size as home range overlap
# not all individuals have home range data
AIm18 <- AIm18[row.names(AIm18) %in% row.names(hrom18), colnames(AIm18) %in% colnames(hrom18)]
AIm19 <- AIm19[row.names(AIm19) %in% row.names(hrom19), colnames(AIm19) %in% colnames(hrom19)]
AIm20 <- AIm20[row.names(AIm20) %in% row.names(hrom20), colnames(AIm20) %in% colnames(hrom20)]
AIm21 <- AIm21[row.names(AIm21) %in% row.names(hrom21), colnames(AIm21) %in% colnames(hrom21)]
AIm22 <- AIm22[row.names(AIm22) %in% row.names(hrom22), colnames(AIm22) %in% colnames(hrom22)]

um18 <- um18[row.names(um18) %in% row.names(hrom18), colnames(um18) %in% colnames(hrom18)]
um19 <- um19[row.names(um19) %in% row.names(hrom19), colnames(um19) %in% colnames(hrom19)]
um20 <- um20[row.names(um20) %in% row.names(hrom20), colnames(um20) %in% colnames(hrom20)]
um21 <- um21[row.names(um21) %in% row.names(hrom21), colnames(um21) %in% colnames(hrom21)]
um22 <- um22[row.names(um22) %in% row.names(hrom22), colnames(um22) %in% colnames(hrom22)]

wm18 <- wm18[row.names(wm18) %in% row.names(hrom18), colnames(wm18) %in% colnames(hrom18)]
wm19 <- wm19[row.names(wm19) %in% row.names(hrom19), colnames(wm19) %in% colnames(hrom19)]
wm20 <- wm20[row.names(wm20) %in% row.names(hrom20), colnames(wm20) %in% colnames(hrom20)]
wm21 <- wm21[row.names(wm21) %in% row.names(hrom21), colnames(wm21) %in% colnames(hrom21)]
wm22 <- wm22[row.names(wm22) %in% row.names(hrom22), colnames(wm22) %in% colnames(hrom22)]

## MRQAP

mrqap18u <- mrqap.dsp(um18 ~ AIm18 + hrom18, randomisations=1000)
mrqap18u
mrqap18w <- mrqap.dsp(wm18 ~ AIm18 + hrom18, randomisations=1000)
mrqap18w 

mrqap19u <- mrqap.dsp(um19 ~ AIm19 + hrom19, randomisations=1000)
mrqap19u # prints results
mrqap19w <- mrqap.dsp(wm19 ~ AIm19 + hrom19, randomisations=1000)
mrqap19w # prints results

mrqap20u <- mrqap.dsp(um20 ~ AIm20 + hrom20, randomisations=1000)
mrqap20u # prints results
mrqap20w <- mrqap.dsp(wm20 ~ AIm20 + hrom20, randomisations=1000)
mrqap20w # prints results

mrqap21u <- mrqap.dsp(um21 ~ AIm21 + hrom21, randomisations=1000)
mrqap21u # prints results
mrqap21w <- mrqap.dsp(wm21 ~ AIm21 + hrom21, randomisations=1000)
mrqap21w # prints results

mrqap22u <- mrqap.dsp(um22 ~ AIm22 + hrom22, randomisations=1000)
mrqap22u # prints results
mrqap22w <- mrqap.dsp(wm22 ~ AIm22 + hrom22, randomisations=1000)
mrqap22w # prints results
