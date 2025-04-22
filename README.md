
# sNuConv

A bulk RNA-seq deconvolution method trained on single-nucleus RNA-seq data.
See the manuscript published in **[iScience](https://www.sciencedirect.com/science/article/pii/S2589004224015931?via%3Dihub)**

**iMPORTANT Erratum statement regarding the published manuscript:**
The authors have discovered an error in identification of sample 3313 sent for snRNA-seq,
such, that it was from a different donor than true-sample 3313 sent for bulk RNA-seq. Hence, in the first cohort, only 6, not 7, hVAT samples were analyzed in parallel by bulk and snRNA-seq.
Re-analysis of the data without sample 3313 showed slightly stronger results than originally reported, thus having no material effect on the results and conclusions of the study

**This Github repo dose not include Sample 3313 anymore**



### Table of Contents    
- [Dependencies](#dependencies)
- [Predicting hVAT\hSAT cell-type Proportions](#Predicting)
  - [Mapping Bulk RNASeq](#Mapping)
  - [Using Docker](#using-docker)
  - [No Conda [Not Recommended]](#no-conda-not-recommended)
- [Tutorial](#tutorial)
- [Contact](#contact)

&nbsp;  
&nbsp;
&nbsp;  


## Dependencies
1. Clone the git repo:
    ```Bash
      bash
      git clone https://github.com/bioinfo-core-BGU/sNuConv.git
      cd sNuConv
    ```  
2. Install mamba and Create the **Scaden** conda environment:
    ```Bash
       conda install conda-forge::mamba
       mamba create -f Scaden_environment.yml
    ```  
3. Create the **RNASeq** conda environment:
   ```Bash
       mamba create -n RNASeq bioconda::rsem bioconda::star
    ```  
   
4. Download the Human Genome and annotation. **It is very important to use this specific versions of genome and annotation since these where used for training Scaden**
    ```Bash
       mkdir Genome
	   cd Genome
	   wget http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz
	   wget http://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	   gzip -d *.gz
	   cd ..
    ```  
## Predicting

  ### Mapping
  Using RSEM to map the Bulk RNASeq Samples to the Human Genome v**GRCh38.100**
  1. Activate **RNASeq** Conda environment:
    ```Bash
        conda activate RNASeq
    ```  
    If it dose not work:
    ```Bash
        source activate RNASeq
    ```  
  2. Prepare reference 
     ```Bash
        rsem-prepare-reference \
            --gtf Genome/Homo_sapiens.GRCh38.100.gtf \
            --star --star-path $CONDA_PREFIX/bin  \
            Genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa  \
            Genome/reference 
     ```  
  3. Map Samples to the reference genome (Do it for all your samples):
     ```Bash
        mkdir Samples
        mkdir Samples/Sample1
        
        rsem-calculate-expression \
        --append-names \
        --estimate-rspd \
        --keep-intermediate-files \
        --output-genome-bam \
        --forward-prob 0.5 \
        -p 10 \
        --star --star-path $CONDA_PREFIX/bin  \
        Sample1.fastq \
        Genome/reference \
        Samples/Sample1 \
        > Samples/Sample1/Sample1_RSEM.out
     ```  
  4. After Mapping all samples, Marge all Samples count data in to one file:
     ```Bash
        rsem-generate-data-matrix  \
           $(find Samples/ -type f -name "*.genes.results")  \
           > Samples/GeneMat.results 
     ```  
  
  ### Predicting
  1. Preprocess:
     ```Bash
        mkdir Predicting_hSAT
		
        Rscript scripts/Preprocess.R \
			--bulk_source RSEM \
			--bulk_counts Samples/GeneMat.results \
			--model_dir models/hSAT/Scaden_model \
			--outDir Predicting_hSAT/ 
			
		mkdir Predicting_hVAT
		
        Rscript scripts/Preprocess.R \
			--bulk_source RSEM \
			--bulk_counts Samples/GeneMat.results \
			--model_dir models/hVAT/Scaden_model \
			--outDir Predicting_hVAT/ 
			
     ```  
  2. Predict
    ```Bash
	   scaden predict \
		--outname Predicting_hSAT//predictions.txt \
		--model_dir models/hSAT/Scaden_model \
		Predicting_hSAT//bulk_data.txt   
		
		scaden predict \
		--outname Predicting_hVAT//predictions.txt \
		--model_dir models/hVAT/Scaden_model \
		Predicting_hVAT//bulk_data.txt   
	 
	```
  3. Correct Predictions
    ```Bash
	   Rscript scripts/CorrectPredictions.R \
		--model models/hSAT/cell_type_regression_model/Model.txt \
		--predictions Predicting_hSAT//predictions.txt \
		--outDir Predicting_hSAT/ 
		
		Rscript scripts/CorrectPredictions.R \
		--model models/hVAT/cell_type_regression_model/Model.txt \
		--predictions Predicting_hVAT//predictions.txt \
		--outDir Predicting_hVAT/ 
	 
	```


## Usage
A sNuConv workflow consists of four major steps:
* Generating per-gene regression model
* Pseudo-bulk simulation
* Deep-learning training using [Scaden](https://github.com/KevinMenden/scaden)
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
