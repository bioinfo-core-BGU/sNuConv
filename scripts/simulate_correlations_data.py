import sys
import pandas as pd
import numpy as np
import datetime

if len(sys.argv) != 7:
    sys.exit("Usage: python simulate_correlations_data.py <snuc_counts.txt> <celltypes.txt> <sample_proportions.txt> <n_cells> <n_simulations> <output_folder>")

started = datetime.datetime.now()
# load the datasets
print("----------------------------------------------------------------")
print("(%s): Reading counts data..." % datetime.datetime.now().strftime("%H:%M:%S"))
data_x = pd.read_table(sys.argv[1], index_col=0, dtype=np.float32)
print("(%s): Done" % datetime.datetime.now().strftime("%H:%M:%S"))
print("----------------------------------------------------------------")
data_y = pd.read_table(sys.argv[2])
sample_props = pd.read_table(sys.argv[3], index_col=0)
sample_size = int(sys.argv[4])
simulations = int(sys.argv[5])
out_folder = sys.argv[6]

# Extract celltypes
celltypes = list(set(data_y["Celltype"].tolist()))
df = []
for sample in sample_props.index:
    print("(%s): Started sample '%s'" % (datetime.datetime.now().strftime("%H:%M:%S"), sample))
    fracs = np.array(sample_props.loc[sample])
    samp_fracs = np.multiply(fracs, sample_size)
    samp_fracs = list(map(int, samp_fracs))
    fracs_df = pd.DataFrame([samp_fracs], index=['Cells'], columns=[sample_props.columns])
    print(fracs_df)
    print("------------------------------------------------------------------------------------------------")
    
    # Simulate for each sample
    sim_x = []
    for j in range(simulations):
        artificial_samples = []
        for i in range(len(celltypes)):
            ct = celltypes[i]
            cells_sub = data_x.loc[np.array(data_y["Celltype"] == ct), :]
            ct_index = np.where(sample_props.columns==ct)[0][0]
            cells_fraction = np.random.randint(0, cells_sub.shape[0], samp_fracs[ct_index])
            cells_sub = cells_sub.iloc[cells_fraction, :]
            artificial_samples.append(cells_sub)
        df_samp = pd.concat(artificial_samples, axis=0)
        df_samp = df_samp.sum(axis=0)
        sim_x.append(df_samp)
    
    simx_df = pd.concat(sim_x, axis=1)
    simx_df.to_csv(out_folder + '/' + sample + '.csv')
    
    # Calculate mean expression for each gene
    df.append(simx_df.mean(axis=1))

# Concat samples and save dataframe
df_mean = pd.concat(df, axis=1)
df_mean.columns = [sample_props.index]
df_norm = (df_mean / df_mean.sum())*1000000
df_norm.to_csv(out_folder + '/correlations_simulation_data.csv')

# Running time
delta = datetime.datetime.now() - started
print("Total running time of the script: %s" % str(delta).split('.')[0])