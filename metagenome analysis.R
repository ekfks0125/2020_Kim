## This code base on DADA2 Pipline tutorial (ver. 1.12)
## https://benjjneb.github.io/dada2/tutorial.html
##                                    the code written by Gyeongjun Cho
##=====================================================================
#===================
# 1. Getting ready
#===================
# DADA2, ggplot2 library loading
library(dada2)
library(DECIPHER)
library(doParallel)
library(ggplot2)

# Input raw files directory
path <- "raw_fastq_files"

# Extract raw files list
list.files(path)

# Foward reads list
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))

# Reverse reads list
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))

# Extraction sample names from file list
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Apply sample name
names(fnFs) <- sample.names 
names(fnRs) <- sample.names

# Sample order
sample_ord <- c("1UT","1F","1FG","1FS4","1FGS","1G","1S4",
                "3UT","3F","3FG","3FS4","3FGS","3G","3S4",
                "7UT",     "7FG","7FS4","7FGS","7G","7S4",
                "10UT",   "10FG","10FS4","10FGS","10G","10S4")


#===================================
# 2. Inspect read quality profiles
#===================================
# Visulaization raw foward qulity
raw_foward_Q <- plotQualityProfile(fnFs[sample_ord])
raw_foward_Q

# Visulaization raw reverse qulity
raw_reverse_Q <- plotQualityProfile(fnRs[sample_ord])
raw_reverse_Q

#=====================
# 3. Filter and trim
#=====================
# Place filtered files in filtered/ subdirectory
filt_path <- "filtered_fastq_files"
filtFs <- file.path(filt_path , paste0(sample.names, "_F_filt.fastq.gz")) # set filtered foward path info
filtRs <- file.path(filt_path , paste0(sample.names, "_R_filt.fastq.gz")) # set filtered reverse path info

# Set sample names for filtered
names(filtFs) <- sample.names 
names(filtRs) <- sample.names  

out <- filterAndTrim(fwd = fnFs, # foward reads
                     filt = filtFs, # filterd foward reads save info
                     rev = fnRs, # reverse reads 
                     filt.rev = filtRs, # filterd reverse reads save info
                     truncLen=c(300,270), # trim to Qulity score 30 or higher
                     trimLeft = c(nchar("GTGYCAGCMGCCGCGGTAA"),
                                  nchar("GGACTACNVGGGTWTCTAAT")), # trim the 515F, 806R primer
                     maxN = 0, maxEE = c(1,1), rm.phix = TRUE, # Other statistical quality check criteria
                     n=1e8, # sampling reads number
                     compress = TRUE, # compress
                     multithread = TRUE) # On Windows set multithread=FALSE

# Visulaization filtered foward qulity
filt_foward_Q <- plotQualityProfile(filtFs[sample_ord])
filt_foward_Q


# Visulaization filtered reverse qulity
filt_reverse_Q <- plotQualityProfile(filtRs[sample_ord])
filt_reverse_Q

#=========================
# 4. Learn the error rates
#=========================
# This step took about 5 hours when 1700X (16 threads) and 3000Mhz DDR4 32Gib ram.
# If you want to save time, reduce the nbase value.
errF <- learnErrors(filtFs, multithread=TRUE, nbases = 1e+10)
errR <- learnErrors(filtRs, multithread=TRUE, nbases = 1e+10)
errF_plot <- plotErrors(errF, nominalQ=TRUE)
errR_plot <- plotErrors(errF, nominalQ=TRUE)
errF_plot
errR_plot

#=====================================================================
# 5. Sample Inference using Divisive Amplicon Denoising Algorithm (DADA)
#=====================================================================
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#======================
# 6. Merge paired reads
#======================
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#=========================
# 7. Construct sequence table
#=========================
seqtab <- makeSequenceTable(mergers)

#====================
# 8. Remove chimeras
#====================
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#===============================
# 9. Assign taxonomy with IDTAXA
#===============================
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet
load("IDTAXA SILVA DB/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
rm(trainingSet) # remove SILVA DB on ram for ram management

#==================================
# 10. none-bacterial OTUs removing
#==================================
taxid_OTUtable <- cbind(data.frame(taxid), as.data.frame(t(seqtab.nochim)[,sample_ord]))
taxid_OTUtable <- subset(taxid_OTUtable, domain == "Bacteria")
taxid_OTUtable <- subset(taxid_OTUtable, order != "Chloroplast")
taxid_OTUtable <- subset(taxid_OTUtable, family != "Mitochondria")

#======================================
# 11. Track reads through the pipeline
#======================================
getN <- function(x) sum(getUniques(x))

inputF <- raw_foward_Q$layers[[6]]$data$rc
names(inputF) <- raw_foward_Q$layers[[6]]$data$label

inputR <- raw_reverse_Q$layers[[6]]$data$rc
names(inputR) <- raw_reverse_Q$layers[[6]]$data$label

filteredF <- filt_foward_Q$layers[[6]]$data$rc
names(filteredF) <- filt_foward_Q$layers[[6]]$data$label

filteredR <- filt_reverse_Q$layers[[6]]$data$rc
names(filteredR) <- filt_reverse_Q$layers[[6]]$data$label

inoutput <- data.frame(inputF,inputR,filteredF,filteredR)

track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
track <- track[sample_ord,]
tax_Bacteria <- colSums(taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
track <- cbind(inputF, inputR, filteredF, filteredR, track, tax_Bacteria)
rownames(track) <- sample_ord
write.csv(track, "Reads_track.csv")

#================
# 12. output files
#================
registerDoParallel(cores=16) # CPU threads to multithreads

# OTUs ID generation
OTUsID <- foreach(i = 1:dim(taxid_OTUtable)[1], .combine="rbind")%dopar%{
  OTUsN <- paste("OTU", i, sep="_")
  print(OTUsN)
}
OTUsID <- as.vector(OTUsID)  
OTUsID

# OTU fasta files
fasta <- foreach(i= 1:dim(taxid_OTUtable)[1], .combine="rbind")%dopar%{
  match <- rbind(OTUsID[i], rownames(taxid_OTUtable)[i])
  print(match)
} 
fasta <- as.vector(fasta)
fasta <- as.data.frame(gsub("^OTU_", ">OTU_",fasta))
colnames(fasta) <- NA
write.table(fasta, file = "OTUs_sequence.fasta", col.names=FALSE, row.names = FALSE, quote = FALSE)

# OTU fasta based phylogenetic tree
#OTU_seq <- c(`OTU seq` = "basic_info_and_analysis_results/1_Microbiome_analysis_step1_results/OTUs_sequence.fasta")
#dbConn <- dbConnect(SQLite(), ":memory:")
#Seqs2DB(OTU_seq, type="FASTA", dbFile = dbConn,"")

#x <- dbGetQuery(dbConn,   "select description from Seqs")$description
#ns <- unlist(lapply(strsplit(x,      split=" "),   FUN=`[`,   1L))
#Add2DB(myData=data.frame(identifier=ns,      stringsAsFactors=FALSE),   dbFile=dbConn)
#dna <- SearchDB(dbConn,   nameBy="identifier")
#align_dna <- AlignSeqs(dna)

#  #calculate a maximum likelihood tree
#  d <- DistanceMatrix(align_dna, correction="Jukes-Cantor")
#  dend <- IdClusters(d,   method="ML",   type="dendrogram",   myXStringSet=align_dna, processors = 16)

#WriteDendrogram(dend, file = "basic_info_and_analysis_results/1_Microbiome_analysis_step1_results/OTUs_tree.newick")
#dbDisconnect(dbConn)

# OTU counts table
OTU_counts_table <- cbind(OTUsID, taxid_OTUtable[,8:dim(taxid_OTUtable)[2]])
row.names(OTU_counts_table) <- NULL
colnames(OTU_counts_table)[1] <- "#OTU ID"
write.table(OTU_counts_table, file = "OTU_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# OTU tax table
OTU_tax_table <- cbind(OTUsID, taxid_OTUtable[,1:7])
row.names(OTU_tax_table) <- NULL
colnames(OTU_tax_table)[1] <- "#OTU ID"
write.table(OTU_tax_table, file = "OTU_tax_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# sequence tax counts table
OTU_seq_tax_counts_table <- cbind(as.vector(row.names(taxid_OTUtable)), OTUsID, taxid_OTUtable)
row.names(OTU_seq_tax_counts_table) <- NULL
write.table(OTU_seq_tax_counts_table, file = "OTU_seq_tax_counts_table.tsv", col.names=TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

#save env.image
save.image(file = "metagenome.RData")


library(phyloseq)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
OTU <- otu_table(
t(
  read.csv("OTU_counts_table.tsv", sep = "\t",row.names = 1)
),
taxa_are_rows = FALSE
)
row.names(OTU) <- gsub("X","",row.names(OTU))

TAX <- tax_table(
  as.matrix(
    read.csv("OTU_tax_table.tsv", sep = "\t",row.names = 1)
  )
)
# Sample inforamtion
SAM <- sample_data(
  read.csv("dr_tomato_pot_test.csv", row.names = 1)
)
SAM
# To make phyloseq obj
# reads base obj
phy_obj_reads <- phyloseq(OTU,TAX,SAM)#,TRE)
phy_obj_reads

# to modified relative OTU phyloseq obj
phy_obj_relative <- phy_obj_reads
phy_obj_relative@otu_table <- phy_obj_reads@otu_table/rowSums(phy_obj_reads@otu_table)
write.csv(t(phy_obj_relative@otu_table[,paste0("OTU_",1:ncol(phy_obj_relative@otu_table))]), "relative_abundance_by_OTU.csv")


relative_abundance <- plot_bar(phy_obj_relative)$data
relative_abundance$genus <- as.character(relative_abundance$genus)
relative_abundance$genus[is.na(relative_abundance$genus)] <- "not assigned"

relative_abundance <- relative_abundance[as.character(1:dim(relative_abundance)[1]),]
relative_abundance$Sampling.time <- factor(relative_abundance$Sampling.time, levels = c("1 sample", "3 sample", "7 sample", "10 sample"))
relative_abundance$Sample <- factor(relative_abundance$Sample, levels = sample_ord)
relative_abundance$Treatment <- factor(relative_abundance$Treatment, levels = unique(as.character(SAM$Treatment)))

phylum_relative_abundance <- aggregate(Abundance ~ phylum + Sample + Treatment+ Sampling.time, relative_abundance, sum)
phylum_relative_abundance$Sample <- factor(phylum_relative_abundance$Sample, levels = rownames(SAM))

class_relative_abundance <- aggregate(Abundance ~ class + Sample + Treatment+ Sampling.time, relative_abundance, sum)
class_relative_abundance$Sample <- factor(class_relative_abundance$Sample, levels = rownames(SAM))

order_relative_abundance <- aggregate(Abundance ~ order + Sample + Treatment+ Sampling.time, relative_abundance, sum)
order_relative_abundance$Sample <- factor(order_relative_abundance$Sample, levels = rownames(SAM))

family_relative_abundance <- aggregate(Abundance ~ family + Sample + Treatment+ Sampling.time, relative_abundance, sum)
family_relative_abundance$Sample <- factor(family_relative_abundance$Sample, levels = rownames(SAM))

genus_relative_abundance <- aggregate(Abundance ~ genus + Sample + Treatment+ Sampling.time, relative_abundance, sum)
genus_relative_abundance$Sample <- factor(genus_relative_abundance$Sample, levels = rownames(SAM))

# top 10 sorting
# phylum_level
sum_phylum <- aggregate(Abundance ~ phylum, relative_abundance, sum)
sum_phylum <- sum_phylum %>% 
  arrange(desc(Abundance))
top_phylum <- as.vector(sum_phylum[1:10,1])

# class_level
sum_class <- aggregate(Abundance ~ class, relative_abundance, sum)
sum_class <- sum_class %>% 
  arrange(desc(Abundance))
top_class <- as.vector(sum_class[1:10,1])

# order_level
sum_order <- aggregate(Abundance ~ order, relative_abundance, sum)
sum_order <- sum_order %>% 
  arrange(desc(Abundance))
top_order <- as.vector(sum_order[1:10,1])

# family_level
sum_family <- aggregate(Abundance ~ family, relative_abundance, sum)
sum_family <- sum_family %>% 
  arrange(desc(Abundance))
top_family <- as.vector(sum_family[1:10,1])

# genus_level
sum_genus <- aggregate(Abundance ~ genus, relative_abundance, sum)
sum_genus <- sum_genus %>% 
  arrange(desc(Abundance))
top_genus <- as.vector(sum_genus[1:10,1])


# taxa top bar graph
top_bar <- function(df, fill_label="Top taxa", top_toxa_names){
  levels(df$Cycling) <- levels(df$Cycling)[c(1,3:10,2)]
  
  df$top_taxa <- df[,1]
  df$top_taxa <- factor(df$top_taxa, levels = c(top_toxa_names,"etc."))
  df$top_taxa[is.na(df$top_taxa)] <- "etc."
  
  agg_df <- aggregate(Abundance ~ top_taxa + Sample + Treatment + Sampling.time, df, sum)
  
  p <- ggplot(data = agg_df, mapping = aes(x = Sample, y=Abundance, fill=top_taxa))+
    geom_bar(stat = "identity")+
    scale_fill_manual(values = c(brewer.pal(n=10, name="Paired"),"#999999"))+
    theme_light()+
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
    theme(axis.text.x = element_text(angle = 90),
          panel.grid.major = element_blank())+
    scale_y_continuous(labels = scales::percent, limits = c(0,1.05), expand = c(0,0))+ # fill option
    ylab("Relative abundance")+
    xlab("")+
    labs(fill=fill_label)
}

top_10_phylum_bar <- top_bar(phylum_relative_abundance, fill_label = "Top 10 phylum", top_toxa_names = top_phylum )
top_10_phylum_bar + facet_grid(.~Sampling.time, scale="free", space='free_x')
top_10_phylum_bar + facet_grid(.~Treatment, scale="free", space='free_x')

top_10_class_bar <- top_bar(class_relative_abundance, fill_label = "Top 10 class", top_toxa_names = top_class )
top_10_class_bar + facet_grid(.~Sampling.time, scale="free", space='free_x')
top_10_class_bar + facet_grid(.~Treatment, scale="free", space='free_x')

top_10_order_bar <- top_bar(order_relative_abundance,  fill_label = "Top 10 order", top_toxa_names = top_order ) 
top_10_order_bar + facet_grid(.~Sampling.time, scale="free", space='free_x')
top_10_order_bar + facet_grid(.~Treatment, scale="free", space='free_x')

top_10_family_bar <- top_bar(family_relative_abundance, fill_label = "Top 10 family", top_toxa_names = top_family )
top_10_family_bar + facet_grid(.~Sampling.time, scale="free", space='free_x')
top_10_family_bar + facet_grid(.~Treatment, scale="free", space='free_x')

top_10_genus_bar <- top_bar(genus_relative_abundance, fill_label = "Top 10 genus", top_toxa_names = top_genus )
top_10_genus_bar+ facet_grid(.~Sampling.time, scale="free", space='free_x')
top_10_genus_bar + facet_grid(.~Treatment, scale="free", space='free_x')
