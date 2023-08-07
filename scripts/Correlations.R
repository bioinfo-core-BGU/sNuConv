library(optparse)
library(dplyr)
library(Biobase)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--bulk_counts"), type="character", default = NA,
              help="Path to Bulk-RNA expression table", metavar = "character"),
  make_option(c("--pseudo_counts"), type="character", default = NA,
              help="Path to Pseudo Bulk-RNA counts CSV file", metavar = "character"),
  make_option(c("--method"), type="character", default = 'pearson',
              help="which correlation coefficient (or covariance) is to be computed (Default is pearson)", metavar = "character"),
  make_option(c("--corr_cutoff"), type="numeric", default = 0.6,
              help="Only use genes with correlations above the cutoff (Default is 0.6)", metavar = "character"),
  make_option(c("--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,'logs.txt',sep='')
cat(paste('[',Sys.time(),']: Correlations',sep=''), file=logfile, sep='\n')
writeLog <- function(msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$bulk_counts)) {
  if (file.exists(opt$bulk_counts)) {
    writeLog(paste("Importing Bulk counts: ",opt$bulk_counts,sep=''))
    bulk_counts = as.matrix(read.table(opt$bulk_counts, row.names = 1, header = TRUE, sep = '\t'))
  } else {
    writeLog(paste("ERROR: Bulk-RNA counts file does not exists."))
    stop()
  }
} else {
  writeLog(paste("ERROR: Bulk-RNA counts file must be specified [--bulk_counts]"))
  stop()
}
if (!is.na(opt$pseudo_counts)) {
  if (file.exists(opt$pseudo_counts)) {
    writeLog(paste("Importing Pseudo counts: ",opt$pseudo_counts,sep=''))
    pseudo_counts = as.matrix(read.csv(opt$pseudo_counts, row.names = 1, header = TRUE))
  } else {
    writeLog(paste("ERROR: Pseudo Bulk-RNA counts file does not exists."))
    stop()
  }
} else {
  writeLog(paste("ERROR: Pseudo Bulk-RNA counts CSV file must be specified [--pseudo_counts]"))
  stop()
}
if (opt$method %in% c('pearson','kendall','spearman')) {
  writeLog(paste("Correlation Coefficient: ",opt$method,sep=''))
} else {
  writeLog(paste("ERROR: Unknown method '",opt$method,"', Please use pearson/kendall/spearman",sep=''))
  stop()
}
if (is.na(opt$outDir)) {
  writeLog(paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}

# Prepare data
common_genes <- intersect(rownames(bulk_counts),rownames(pseudo_counts))
writeLog(paste("Found ",length(common_genes)," shared genes between Bulk & Pseudo datasets",sep=''))
bulk_counts = as.matrix(bulk_counts[common_genes,order(colnames(bulk_counts))])
pseudo_counts = as.matrix(pseudo_counts[common_genes,order(colnames(pseudo_counts))])
if (!(all(colnames(bulk_counts)==colnames(pseudo_counts)))) {
  writeLog(paste("ERROR: The sample names are different in Bulk & Pseudo counts datasets"))
  stop()
}

# Correlations
writeLog(paste("Calculating correlations..."))
df_corr = data.frame()
for (gene in rownames(bulk_counts)) {
  true_vec = bulk_counts[gene,] %>% as.vector()
  pseudo_vec = pseudo_counts[gene,] %>% as.vector()
  gene_cor = cor(true_vec, pseudo_vec, method = opt$method)
  gene_lm = lm(true_vec ~ pseudo_vec)
  new_gene_df = data.frame('Gene'=gene, 'Correlation'=gene_cor, 'abs_cor'=abs(gene_cor), 
                           'Slope'=gene_lm$coefficients[2], 'Intercept'=gene_lm$coefficients[1])
  df_corr = rbind(df_corr, new_gene_df)
}
rownames(df_corr) = df_corr[,'Gene']
df_corr <- df_corr[!is.na(df_corr$abs_cor),]
genes_to_correct = df_corr$Gene[df_corr$abs_cor>=opt$corr_cutoff]
writeLog(paste("Found ",length(genes_to_correct)," genes with absolute correlation above ",opt$corr_cutoff,sep=''))
df_corr = df_corr[genes_to_correct,]

# End of run
writeLog(paste("Saving results..."))
write.table(df_corr, paste(opt$outDir,'df_corr',sep=''), sep = '\t', quote = FALSE, row.names = FALSE)
writeLog(paste("Finished"))