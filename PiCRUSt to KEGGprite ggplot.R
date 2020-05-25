### Reproducing the categorize by function (level 3) functionality in plain-text tables.
### Doing this because adding a column of KEGG Pathways to a table and then converting
### that table to BIOM is difficult.

categorize_by_function_l3 <- function(in_ko, kegg_brite_mapping) {
  # Function to create identical output as categorize_by_function.py script,
  # but with R objects instead of BIOM objects in Python.
  # Input KO table is assumed to have rownames as KOs and sample names as columns.
  
  out_pathway <- data.frame(matrix(NA, nrow=0, ncol=(ncol(in_ko) + 1)))
  
  colnames(out_pathway) <- c("pathway", colnames(in_ko))
  
  for(ko in rownames(in_ko)) {
    
    # Skip KO if not in KEGG BRITE mapping df
    # (this occurs with newer KOs that weren't present in PICRUSt1).
    if(! ko %in% rownames(kegg_brite_mapping)) {
      next
    }
    
    pathway_list <- strsplit(kegg_brite_mapping[ko, "metadata_KEGG_Pathways"], "\\|")[[1]]
    
    for(pathway in pathway_list) {
      
      pathway <- strsplit(pathway, ";")[[1]][3]
      
      new_row <- data.frame(matrix(c(NA, as.numeric(in_ko[ko,])), nrow=1, ncol=ncol(out_pathway)))
      colnames(new_row) <- colnames(out_pathway)
      new_row$pathway <- pathway
      out_pathway <- rbind(out_pathway, new_row)
    }
    
  }
  
  out_pathway = data.frame(aggregate(. ~ pathway, data = out_pathway, FUN=sum))
  
  rownames(out_pathway) <- out_pathway$pathway
  
  out_pathway <- out_pathway[, -which(colnames(out_pathway) == "pathway")]
  
  if(length(which(rowSums(out_pathway) == 0)) > 0) {
    out_pathway <- out_pathway[-which(rowSums(out_pathway) == 0), ]
  }
  
  return(out_pathway)
  
}
#=============================================================================================================================================
kegg_brite_map <- read.table("/home/heuklang/다운로드/picrust1_KO_BRITE_map.tsv", header = T, sep="\t", quote= "", stringsAsFactors = F, comment.char="", row.names=1)

reading_ko <- read.table("/home/heuklang/prebiotic/picrust2/KO_metagenome_out/pred_metagenome_unstrat.tsv",header=TRUE, sep="\t", row.names=1) 
reading_ko

ko_L3 <- categorize_by_function_l3(reading_ko, kegg_brite_map)

colnames(ko_L3) <- gsub("^X","", colnames(ko_L3) )
colnames(ko_L3)

ko_L3

ko_L3_integer<- ko_L3

for(i in 1:length(colnames(ko_L3_integer))){
  ko_L3_integer[, i] <- as.integer(ko_L3_integer[, colnames(ko_L3_integer)[i]])
}


for_deseq <- data.frame(row.names(ko_L3))
colnames(for_deseq) <- "funtion"
for_deseq <- cbind(for_deseq, ko_L3_integer, row.names=NULL)
for_deseq


metadata <- read.csv("/home/heuklang/prebiotic/meta.tsv", sep="\t")
row.names(metadata) <- metadata[,1]
row.names(metadata)
require(DESeq2)
for_deseq
dds <- DESeqDataSetFromMatrix(countData=for_deseq,
                              colData = metadata,
                              design = ~TG, tidy=T)


dds <- DESeq(dds, test="LRT", reduced= ~ 1)

years_Asp_vs_Glu <- results(dds, c("TG", "Asp","Glu"))$log2FoldChange
years_Ck_vs_Glu <- results(dds, c("TG", "Ck","Glu"))$log2FoldChange
years_Asp_vs_ck <- results(dds, c("TG", "Asp","Ck"))$log2FoldChange
LRT_model_results <- cbind(results(dds)[,-2], years_Asp_vs_Glu, years_Ck_vs_Glu, years_Asp_vs_ck)

colnames(LRT_model_results[,6:8]) <- c("Log2FoldChange_AspvsGlu", "Log2FoldChange_CkvsGlu", "Log2FoldChange_AspvsCK")

write.csv(LRT_model_results, "/home/heuklang/prebiotic/PiCRUST KEGG brite, DESeq LRT analysis.csv")

cutoff_padj <- 0.001
signi_LRT_model_results <- subset(LRT_model_results, padj < cutoff_padj)
signi_LRT_model_results


PiCRUST2_KEGG_brite_of_DESeq_LRT_signi<- ko_L3[rownames(signi_LRT_model_results),]
PiCRUST2_KEGG_brite_of_DESeq_LRT_signi$function_group <- rownames(PiCRUST2_KEGG_brite_of_DESeq_LRT_signi)
PiCRUST2_KEGG_brite_of_DESeq_LRT_signi$padj <- signi_LRT_model_results$padj

library(reshape2)
for_ggplot_signi_KEGG <- melt(PiCRUST2_KEGG_brite_of_DESeq_LRT_signi, id=c("function_group","padj"))
for_ggplot_signi_KEGG 
colnames(for_ggplot_signi_KEGG) <- c("function_group", "padj","sample", "PiCRUST2_predicted_KEGG_brite") 
for_ggplot_signi_KEGG$sample
for_ggplot_signi_KEGG$amino <- gsub("^[0-9Sampling]{1,}_","",for_ggplot_signi_KEGG$sample)
for_ggplot_signi_KEGG$amino <- gsub("[0-9].fastq","",for_ggplot_signi_KEGG$amino)
for_ggplot_signi_KEGG$amino

for_ggplot_signi_KEGG

changeSciNot <- function(n){
  a <-c()
  for(i in 1:length(n)){
    output <- format(n, scientific = TRUE, digits = 3) #Transforms the number into scientific notation even if small
    output <- sub("e", "%*%10^", output) #Replace e with 10^
    output <- sub("\\+0?", "", output) #Remove + symbol and leading zeros on expoent, if > 1
    output <- sub("-0?", "-", output) #Leaves - symbol but removes leading zeros on expoent, if < 1
    a <- c(output)
  }
  print(a)
}

for_ggplot_signi_KEGG$Padj_form <- changeSciNot(for_ggplot_signi_KEGG$padj)
for_ggplot_signi_KEGG$Padj_labeller <- paste("italic(p)[adj] ==",for_ggplot_signi_KEGG$Padj_form)
for_ggplot_signi_KEGG$Padj_labeller


require(ggplot2)
P <- ggplot(data = for_ggplot_signi_KEGG,
       mapping = aes(x = amino, y=PiCRUST2_predicted_KEGG_brite, color=amino, fill=amino, group=amino))+
  stat_summary(geom="bar", fun.data=mean_sdl, width=.8, fill="white", position = position_dodge(1.2))+
  stat_summary(geom="errorbar", fun.data=mean_sdl, width=0.4,position = position_dodge(1.2))+
  geom_jitter(stat="identity")+
  facet_wrap(~function_group, scales="free")+
  theme_minimal()+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=1)+
  ggtitle(expression(paste("DESeq Likehood Ratio Test, ", italic(p)[adj] < 0.01)))+
  ylab("KEGG brite categorized PiCRUSt2")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank())+
  geom_text(aes(label = Padj_labeller, x = Inf, y=Inf), hjust=1.2, vjust=2.2 ,parse=T, color="black", size=2.5, alpha=0.15)

P
