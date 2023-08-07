library(optparse)
library(dplyr)
options(digits=10)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--truep"), type="character", default = NA,
              help="Path to sNuc-Seq true proportions txt file", metavar = "character"),
  make_option(c("--predictions"), type="character", default = NA,
              help="Path to predictions output file from Scaden", metavar = "character"),
  make_option(c("--Celltypes"), action="store_true", default = FALSE,
              help="Build model based on celltypes", metavar = "character"),
  make_option(c("--Samples"), action="store_true", default = FALSE,
              help="Build model based on samples", metavar = "character"),
  make_option(c("-o", "--outDir"), type="character", default = NA,
              help="Path to the output directory", metavar = "character")
);

#Get user information
opt_parser = optparse::OptionParser(usage = "usage: %prog [options]", 
                                    option_list=option_list,
                                    epilogue="\n\nAuthor:Gil Sorek");
opt = optparse::parse_args(opt_parser);

# Set log framework
logfile = paste(opt$outDir,'log.txt',sep='')
cat(paste('[',Sys.time(),']: Create Model',sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$truep)) {
  writeLog(logfile, paste("Importing sNuc-Seq true proportions...",sep=''))
  truep <- as.matrix(read.delim(opt$truep,header = TRUE, row.names=1))
} else {
  writeLog(logfile, paste("ERROR: sNuc-Seq true proportions file must be specified [--truep]"))
  stop()
}
if (!is.na(opt$predictions)) {
  writeLog(logfile, paste("Importing predictions file...",sep=''))
  preds <- as.matrix(read.delim(opt$predictions,header = TRUE, row.names=1))
} else {
  writeLog(logfile, paste("ERROR: Scaden's predictions file must be specified [--predictions]"))
  stop()
}
if (!opt$Celltypes && !opt$Samples) {
  writeLog(logfile, paste("ERROR: Model Creation method must be specified [--Celltypes or --Samples]"))
  stop()
}

# Filter and sort data
writeLog(logfile, paste("Filtering and Sorting data...",sep=''))
samples = intersect(rownames(preds),rownames(truep))
cts = intersect(colnames(preds),colnames(truep))
truep <- truep[samples,cts]
preds <- preds[samples,cts]
truep <- truep[order(rownames(truep)),]
truep <- truep[,order(colnames(truep))]
preds <- preds[order(rownames(preds)),]
preds <- preds[,order(colnames(preds))]
samples = rownames(truep)
cts = colnames(truep)
writeLog(logfile, paste("Samples detected: ",paste(samples,collapse=','),sep=''))
writeLog(logfile, paste("Celltypes detected: ",paste(cts,collapse=','),sep=''))

# Build Linear Models
if (opt$Celltypes) {
  writeLog(logfile, paste("Building model using Celltypes...",sep=''))
  lms = data.frame()
  for (ct in cts) {
    true = truep[,ct] %>% as.vector()
    pred = preds[,ct] %>% as.vector()
    lm = lm(true ~ pred)
    df = data.frame('Celltype'=ct, 'Slope'=lm$coefficients[2], 'Intercept'=lm$coefficients[1], 'Rsquare'=summary(lm)$r.squared)
    lms = rbind(lms, df)
  }
} else if (opt$Samples) {
  writeLog(logfile, paste("Building model using Samples...",sep=''))
  lms = data.frame()
  for (sample in samples) {
    true = truep[sample,] %>% as.vector()
    pred = preds[sample,] %>% as.vector()
    lm = lm(true ~ pred)
    df = data.frame('Sample'=sample, 'Slope'=lm$coefficients[2], 'Intercept'=lm$coefficients[1])
    lms = rbind(lms, df)
  }
}

# Save results
writeLog(logfile, paste("Saving results...",sep=''))
write.table(lms, paste(opt$outDir,'Model.txt',sep=''), row.names=FALSE, sep='\t')
writeLog(logfile, paste("Finished",sep=''))
