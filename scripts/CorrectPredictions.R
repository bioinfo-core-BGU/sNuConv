library(optparse)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("--model"), type="character", default = NA,
              help="Path to Model txt file", metavar = "character"),
  make_option(c("--predictions"), type="character", default = NA,
              help="Path to predictions output file from Scaden", metavar = "character"),
  make_option(c("--R_cutoff"), type="numeric", default = 0.8,
              help="Cutoff for correcting cell-type proportions (Default is 0.8)", metavar = "character"),
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
cat(paste('[',Sys.time(),']: Correct Predictions',sep=''), file=logfile, sep='\n')
writeLog <- function(logfile, msg) { cat(paste('(',strftime(Sys.time(), format = "%H:%M:%S"),'): ',msg,sep=''),file=logfile,append=TRUE,sep='\n') }
saveRDS(opt, paste(opt$outDir,'opt.rds',sep=''))

print("Input var:",quote = F)
print.data.frame(as.data.frame(x = unlist(opt),row.names = names(opt)),right = F,quote = F)

# Import data
if (!is.na(opt$model)) {
  writeLog(logfile, paste("Importing Model...",sep=''))
  model <- read.delim(opt$model,header = TRUE)
  model_source = colnames(model)[1]
  if (!(model_source %in% c('Celltype','Sample'))) {
    writeLog(logfile, paste("ERROR: The model (1st column name) should be based on either Celltype or Sample"))
	stop()
  }
  writeLog(logfile, paste("Correction will be based on: ",model_source,sep=''))
} else {
  writeLog(logfile, paste("ERROR: Model txt file must be specified [--model]"))
  stop()
}
if (!is.na(opt$predictions)) {
  writeLog(logfile, paste("Importing predictions file...",sep=''))
  preds <- read.delim(opt$predictions,header = TRUE, row.names=1)
} else {
  writeLog(logfile, paste("ERROR: Scaden's predictions file must be specified [--predictions]"))
  stop()
}

# Validating model and predictions
writeLog(logfile, paste("Validating model and predictions..."))
if (model_source == 'Celltype') {
  celltypes = model$Celltype
  celltypes = celltypes[order(celltypes)]
  if (!(all(colnames(preds) %in% celltypes))) {
    writeLog(logfile, paste("ERROR: No model found for the following celltypes: ",paste(colnames(preds)[!(colnames(preds) %in% celltypes)],collapse=','),sep=''))
	stop()
  }
  preds = preds[,celltypes]
  samples = rownames(preds)
  if (length(samples)>1) {
    samples = samples[order(samples)]
	preds = preds[samples,]
  }
} else if (model_source == 'Sample') {
  avg_slope = mean(model$Slope)
  avg_intercept = mean(model$Intercept)
}

# Correct Proportions
writeLog(logfile, paste("Correcting Proportions..."))
preds_corrected = data.frame(row.names=samples)
cts_adjusted = c()
if (model_source == 'Celltype') {
  for (ct in celltypes) {
    Rsquare = model$Rsquare[model$Celltype==ct]
	if (Rsquare >= opt$R_cutoff) {
	  cts_adjusted = c(cts_adjusted, ct)
	  slope = model$Slope[model$Celltype==ct]
	  intercept = model$Intercept[model$Celltype==ct]
	  df = data.frame(preds[,ct] * slope + intercept)
	} else {
	  df = data.frame(preds[,ct])
	}
	rownames(df) = samples
	colnames(df) = ct
	preds_corrected = cbind(preds_corrected, df)
  }
} else if (model_source == 'Sample') {
    preds_corrected = preds * avg_slope + avg_intercept
}

# Normalize proportions to 1
preds_corrected[preds_corrected<0] = 0
for (nrow in 1:nrow(preds_corrected)) {
  props = preds_corrected[nrow,] %>% as.numeric
  names(props) = colnames(preds_corrected)
  lock_props = sum(props[cts_adjusted])
  unlocked_props = props[!(names(props) %in% cts_adjusted)]
  unlocked_props = (props[!(names(props) %in% cts_adjusted)]/sum(unlocked_props))*(1-lock_props)
  for (ct in names(props)) {
    if (!(ct %in% cts_adjusted))
	  props[ct] = unlocked_props[ct]
  }
  preds_corrected[nrow,] = props
}


# Save results
writeLog(logfile, paste("Saving results...",sep=''))
write.table(preds_corrected, paste(opt$outDir,'corrected_predictions.txt',sep=''), row.names=TRUE, quote=FALSE, sep='\t')
system(paste("sed -i '1s/^/\t/' ", paste(opt$outDir,'corrected_predictions.txt',sep=''),sep = '') , intern = F)
writeLog(logfile, paste("Finished",sep=''))
