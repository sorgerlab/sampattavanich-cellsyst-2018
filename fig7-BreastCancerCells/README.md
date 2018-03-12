The scripts in this folder reproduce various figure panels related to Figure 7 and S7.
The scripts assume that a number of `.mat` files are available in the rawdata/fixedcell
folder, with one `.mat` file for each cell line and readout.

# Plotting median-vs-iqr FoxO3 translocation 

Run `analysis_median_iqr_rotation.m` to produce Figure 7B.

# Plotting the effect of MEK/AKT inhibition on AKT/ERK

Run `plot_inhib_effect.m` to produce the following figure panels: Figure 7C, S7A, S7B

# DBN learning of ERK/AKT/Foxo3A connectivity

Run `run_sensitivity.m` to produce the following figure panels: Figure 7D, S7E, S7F.
