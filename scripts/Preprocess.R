library(optparse)
library(dplyr)
library(Biobase)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--bulk_counts"), type="character", default = NA,
              help="Path to Bulk-RNA ExpressionSet", metavar = "character"),
  make_option(c("--bulk_source"), type="character", default = NA,
              help="How crated the bulk counts file [options: RSEM,other]", metavar = "character"),
  make_option(c("--model_dir"), type="character", default = NA,
              help="Path to Scaden's model directory", metavar = "character"),
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
cat(paste('[',Sys.time(),']: Preprocess',sep=''), file=logfile, sep='\n')
writeLog <- function(msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$bulk_counts)) {
  writeLog(paste("Importing bulk data..."))
  bulk_counts = read.table(opt$bulk_counts, header = T, sep = '\t', row.names=1)
  if (opt$bulk_source=='RSEM'){
    genes = gsub(".*_","",rownames(bulk_counts))
    bulk_counts = bulk_counts[!(genes %in% genes[duplicated(genes)]),]
    rownames(bulk_counts) = gsub(".*_","",rownames(bulk_counts)) 
  }
  writeLog(paste("Identified ",ncol(bulk_counts)," samples: ",paste(colnames(bulk_counts),collapse=','),sep=''))
} else {
  writeLog(paste("ERROR: bulk counts data must be specified [--bulk_counts]"))
  stop()
}
if (is.na(opt$outDir)) {
  writeLog(paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}
if (!is.na(opt$model_dir)) {
  writeLog(paste("Importing genes list..."))
  genes_file = paste(opt$model_dir,'m256/genes.txt',sep='')
  if (file.exists(genes_file)) {
    genes = read.table(genes_file, header = T, row.names=1)
	genes = genes[,'X0']
	writeLog(paste("Found ",length(genes)," genes",sep=''))
  } else {
    writeLog(paste("ERROR: Couldn't find genes file in Scaden's model directory"))
	stop()
  }
} else {
  writeLog(paste("ERROR: Scaden's model directory must be specified [--outDir]"))
  stop()
}

# Normalize data
writeLog(paste("Normalizing data..."))
bulk_counts = apply(bulk_counts, MARGIN = 2, FUN=function(x) (x/sum(x))*1000000)

# Filter genes
writeLog(paste("Filtering genes..."))
bulk_counts = bulk_counts[rownames(bulk_counts) %in% genes,]

# End of run
write.table(bulk_counts, file=paste(opt$outDir,'bulk_data.txt',sep=''), sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'bulk_data.txt',sep=''),sep = '') , intern = F)
writeLog(paste("Finished"))
