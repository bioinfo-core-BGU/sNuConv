
# sNuConv

A bulk RNA-seq deconvolution method trained on single-nucleus RNA-seq data.
## Usage
A sNuConv workflow consists of four major steps:
* Generating per-gene regression model
* Pseudo-bulk simulation
* Deep-learning training using Scaden
* Computing cell-type regression model

If training data in ExpressionSet (ESET) is already available, you can start at creating the datasets for Scaden. Otherwise you will first have to process snRNA-seq and bulk RNA-seq datasets to generate a training dataset as ExpressionSet (ESET).
Scaden requires the data to be in a specific format.

First, create the datasets in the required format for Scaden:

```
Rscript CreateMatrices.R \
	--snuc_ident_col orig.ident \
	--snuc_celltype_col Cell_type \
	--bulk_ident_col Sample_name \
	--snuc_eset snuc_eset.rds \
	--bulk_eset bulk_eset.rds \
	--outDir CreateMatrices/
```

This generates the files "snuc_counts.txt", "snuc_celltypes.txt", "snuc_truep.txt" and "bulk_data_training.txt" in the "CreateMatrices" directory. Next, you can run pseudo-bulk simulation using Scaden:
```
scaden simulate \
	-n 1000 \
	-c 1000 \
	-p "training_data" \
	--pattern "*_counts.txt" \
	-d CreateMatrices/ \
	-o Pseudo_Bulk_Simulation/
```
This generates 1000 samples of training data in the "Pseudo_Bulk_Simulation" directory. The file you need is called "training_data.h5ad".

In parallel, we simulate the data to create the per-gene regression model using the true snRNA-seq proportions
```
python simulate_correlations_data.py \
    CreateMatrices/snuc_counts.txt \
	CreateMatrices/snuc_celltypes.txt \
    CreateMatrices/snuc_truep.txt \
    10000 \
	30 \
	Simulate_Correlations_data/
```
This generates the file "correlations_simulation_data.csv" in the "Simulate_Correlations_data" directory.

Next, we compute the per-gene regression model:
```
Rscript Correlations.R \
	--corr_cutoff 0.6 \
	--method spearman \
	--bulk_counts CreateMatrices/bulk_data_training.txt \
	--pseudo_counts Simulate_Correlations_data/correlations_simulation_data.csv \
	--outDir Calculate_Correlations/
```
As a result, you should now have a file called "df_corr", a data frame object by R with the genes regression model passing the 'corr_cutoff' cutoff.

Next, we correct the pseudo-bulk simulation data generated in the "Pseudo_Bulk_Simulation" step.
```
python CorrectCounts.py \
	Calculate_Correlations/df_corr \
	Pseudo_Bulk_Simulation/training_data.h5ad \
	Correct_Pseudo_Counts/
```
You should now have a corrected pseudo-bulk simulation dataset "training_data_corrected.h5ad" ready to be processed by Scaden.
We continue according to Scaden's documentation for data processing and training
```
scaden process \
	--var_cutoff 0 \
	--processed_path Scaden_process/processed.h5ad \
	Correct_Pseudo_Counts/training_data_corrected.h5ad \
	CreateMatrices/bulk_data_training.txt
```
```
scaden train \
	--steps 5000 \
    --model_dir Scaden_train/ \
	Scaden_process/processed.h5ad
```
Finally, you can use the trained model to perform prediction:
```
scaden predict \
	--outname Scaden_predict_training_data/predictions.txt \
	--model_dir Scaden_train/ \
	bulk_data_training.txt
```

As a last step, we compute cell-type regression model using the "snuc_truep.txt" and the prediction results created by Scaden
```
Rscript CreateModel.R \
	--Celltypes \
	--truep CreateMatrices/snuc_truep.txt \
	--predictions Scaden_predict_training_data/predictions.txt \
	--outDir Create_Model/
```
As a result, we should have a "Model.txt" file which includes a per cell-type regression model with the corresponding $R^2$ value.

For the final corrected predictions, you should use the initial predictions with the "Model.txt" file, and choose on a threshold to use for correction based on the $R^2$ (Default: 0.8)
```
Rscript CorrectPredictions.R \
	--model Create_Model/Model.txt \
	--predictions Scaden_predict_training_data/predictions.txt \
	--outDir Correct_Predictions/
```
Now you should have a file called "corrected_predictions.txt" in the "Correct_Predictions" directory, which contains your corrected cell-type proportions.
