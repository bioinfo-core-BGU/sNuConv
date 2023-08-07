library(optparse)
library(dplyr)
library(Biobase)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--snuc_eset"), type="character", default = NA,
              help="Path to Single-Nucleus RNA ExpressionSet", metavar = "character"),
  make_option(c("--bulk_eset"), type="character", default = NA,
              help="Path to Bulk-RNA ExpressionSet", metavar = "character"),
  make_option(c("--snuc_ident_col"), type="character", default = NA,
              help="feature name in the sNuc ESET meta-data of the cell's sample identity", metavar = "character"),
  make_option(c("--snuc_celltype_col"), type="character", default = NA,
              help="feature name in the sNuc ESET meta-data of the cell's cell-type affiliation", metavar = "character"),
  make_option(c("--bulk_ident_col"), type="character", default = NA,
              help="feature name in the Bulk ESET meta-data of the sample's identity", metavar = "character"),
  make_option(c("--test_samples"), type="character", default = NA,
              help="Which samples to use as test data (Default behaviour will train and predict on all samples)", metavar = "character"),
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
cat(paste('[',Sys.time(),']: CreateMatrices',sep=''), file=logfile, sep='\n')
writeLog <- function(msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$snuc_eset)) {
  if (file.exists(opt$snuc_eset)) {
    writeLog(paste("Importing sNuc ESET: ",opt$snuc_eset,sep=''))
    snuc_eset <- readRDS(opt$snuc_eset)
    if (!is.na(opt$snuc_ident_col)) {
      if (opt$snuc_ident_col %in% colnames(snuc_eset@phenoData@data)) {
        snuc_eset[[opt$snuc_ident_col]] = make.names(snuc_eset[[opt$snuc_ident_col]])
        snuc_samples <- unique(snuc_eset[[opt$snuc_ident_col]])
        snuc_samples <- snuc_samples[order(snuc_samples, decreasing = FALSE)]
        writeLog(paste("Loaded ",length(snuc_samples)," Samples: ",paste(snuc_samples,collapse = ','),sep=''))
      } else {
        writeLog(paste("ERROR: Could not find '",opt$snuc_ident_col,"' feature in sNuc meta-data",sep=''))
        stop()
      }
    } else {
      writeLog(paste("ERROR: Cells identity feature name in sNuc ESET must be specified [--snuc_ident_col]"))
      stop()
    }
    if (!is.na(opt$snuc_celltype_col)) {
      if (!(opt$snuc_celltype_col %in% colnames(snuc_eset@phenoData@data))) {
        writeLog(paste("ERROR: Could not find '",opt$snuc_celltype_col,"' feature in sNuc meta-data",sep=''))
        stop()
      }
    } else {
      writeLog(paste("ERROR: Cells cell-type feature name in sNuc ESET must be specified [--snuc_celltype_col]"))
      stop()
    }
  } else {
    writeLog(paste("ERROR: Single-Nucleus RNA ExpressionSet does not exists."))
    stop()
  }
} else {
  writeLog(paste("ERROR: sNuc ESET must be specified [--snuc_eset]"))
  stop()
}
if (!is.na(opt$bulk_eset)) {
  if (file.exists(opt$bulk_eset)) {
    writeLog(paste("Importing Bulk ESET: ",opt$bulk_eset,sep=''))
    bulk_eset <- readRDS(opt$bulk_eset)
    if (!is.na(opt$bulk_ident_col)) {
      if (opt$bulk_ident_col %in% colnames(bulk_eset@phenoData@data)) {
        bulk_eset[[opt$bulk_ident_col]] = make.names(bulk_eset[[opt$bulk_ident_col]])
        bulk_samples <- unique(bulk_eset[[opt$bulk_ident_col]])
        bulk_samples <- bulk_samples[order(bulk_samples, decreasing = FALSE)]
        writeLog(paste("Loaded ",length(bulk_samples)," Samples: ",paste(bulk_samples,collapse = ','),sep=''))
        if (length(bulk_samples) != length(snuc_samples)) {
          writeLog(paste("ERROR: Number of sNuc samples (",length(snuc_samples),") is different from the number of Bulk samples (",length(bulk_samples),")",sep=''))
          stop()
        }
        if (!(all(snuc_samples==bulk_samples))) {
          writeLog(paste("ERROR: The sample names are different in sNuc & Bulk ExpressionSets"))
          stop()
        }
      } else {
        writeLog(paste("ERROR: Could not find '",opt$bulk_ident_col,"' feature in Bulk meta-data",sep=''))
        stop()
      }
    } else {
      writeLog(paste("ERROR: Samples identity feature name in Bulk ESET must be specified [--bulk_ident_col]"))
      stop()
    }
  } else {
    writeLog(paste("ERROR: Bulk-RNA ExpressionSet does not exists."))
    stop()
  }
} else {
  writeLog(paste("ERROR: Bulk-RNA ESET must be specified [--bulk_eset]"))
  stop()
}
if (is.na(opt$outDir)) {
  writeLog(paste("ERROR: Output directory must be specified [--outDir]"))
  stop()
}
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))
if (!is.na(opt$test_samples)) {
  writeLog(paste("Splitting dataset into training data and test data..."))
  test_samples = unlist(stringi::stri_split(str = opt$test_samples,fixed = ','))
  test_samples = test_samples[order(test_samples)]
  if (!(all(test_samples %in% snuc_samples))) {
    writeLog(paste("ERROR: Couldn't find the following samples: ",paste(test_samples[which(!(test_samples %in% snuc_samples))], collapse=','),sep=''))
	stop()
  }
  writeLog(paste("The following samples will be used only as test data: ",paste(test_samples,collapse=','),sep=''))
}

# Export true Cell-type proportions
writeLog(paste("Calculating true cell-type proportions..."))
pData = snuc_eset@phenoData@data[,c(opt$snuc_ident_col, opt$snuc_celltype_col)]
props = pData %>% group_by(opt$snuc_ident_col) %>% table %>% prop.table(margin=1)
truep = as.matrix.data.frame(props)
rownames(truep) = make.names(rownames(props))
colnames(truep) = make.names(colnames(props))
write.table(truep, file=paste(opt$outDir,'snuc_truep.txt',sep=''), sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'snuc_truep.txt',sep=''),sep = '') , intern = F)
if (!is.na(opt$test_samples)) {
  training_truep = truep[rownames(truep)[!(rownames(truep) %in% test_samples)],]
  write.table(training_truep, file=paste(opt$outDir,'training_truep.txt',sep=''), sep='\t')
  system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'training_truep.txt',sep=''),sep = '') , intern = F)
  snuc_eset = snuc_eset[,stringr::str_detect(snuc_eset[[opt$snuc_ident_col]],paste(test_samples,collapse='|'),negate=TRUE)]
} else {
  write.table(truep, file=paste(opt$outDir,'training_truep.txt',sep=''), sep='\t')
  system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'training_truep.txt',sep=''),sep = '') , intern = F)
}

# Export sNuc Celltypes
writeLog(paste("Writing sNuc celltypes matrix..."))
snuc_celltypes <- data.frame('Celltype'=make.names(snuc_eset[[opt$snuc_celltype_col]]))
write.table(snuc_celltypes, file=paste(opt$outDir,'snuc_celltypes.txt',sep=''), sep='\t', row.names = FALSE)

# Export sNuc expression counts
writeLog(paste("Writing sNuc expression matrix..."))
snuc_counts = exprs(snuc_eset)
colnames(snuc_counts) <- 1:length(colnames(snuc_counts))
snuc_counts <- t(snuc_counts)
write.table(snuc_counts, file=paste(opt$outDir,'snuc_counts.txt',sep=''), sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'snuc_counts.txt',sep=''),sep = '') , intern = F)

# Export Bulk expression counts
writeLog(paste("Writing Bulk expression matrix..."))
if (!is.na(opt$test_samples)) {
  bulk_eset_training = bulk_eset[,stringr::str_detect(bulk_eset[[opt$bulk_ident_col]],paste(test_samples,collapse='|'),negate=TRUE)]
  bulk_eset_test = bulk_eset[,stringr::str_detect(bulk_eset[[opt$bulk_ident_col]],paste(test_samples,collapse='|'),negate=FALSE)]
} else {
  bulk_eset_training = bulk_eset
  bulk_eset_test = bulk_eset
}
bulk_counts = exprs(bulk_eset_training)
bulk_counts = apply(bulk_counts, MARGIN = 2, FUN=function(x) (x/sum(x))*1000000)
bulk_counts = bulk_counts[apply(bulk_counts!=0, MARGIN = 1, FUN=any),]
write.table(bulk_counts, file=paste(opt$outDir,'bulk_data_training.txt',sep=''), sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'bulk_data_training.txt',sep=''),sep = '') , intern = F)
test_counts = exprs(bulk_eset_test)
test_counts = apply(test_counts, MARGIN = 2, FUN=function(x) (x/sum(x))*1000000)
write.table(test_counts, file=paste(opt$outDir,'bulk_data_test.txt',sep=''), sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'bulk_data_test.txt',sep=''),sep = '') , intern = F)

# End of run
writeLog(paste("Finished"))