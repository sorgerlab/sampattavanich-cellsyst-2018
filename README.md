# Sampattavanich et al., Cell Systems (2018)

**Source code and example plots for the article "Encoding growth factor identity in the temporal 
dynamics of FoxO3 under combinatorial regulation by ERK and Akt" by Sampattavanich et al., Cell Systems (2018)**

Each piece of source code is provided in a folder containing matlab or python
scripts and all related functions. The main function to generate related
figures are highlighted in **bold**. 

To run this code, users must download related raw data from https://doi.org/10.6084/m9.figshare.c.4026994
and put these in the \\rawdata folder placed at the top-most level of this git
repository folder. users must also install functional data analysis MATLAB package from 
http://www.psych.mcgill.ca/misc/fda/downloads/FDAfuns/Matlab/.

| Figures             | Related files                                                                                                                                                        |
|---------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Fig.1B              | Fig1B/`createParentalVSReporter.m` (Main script)                                                                                                                     |
|                     | Rawdata/parentalVSReporter/\*.\* (download)                                                                                                                          |
| Fig.S1B and S1C     | FigS1BC/`comparedataat15min.m` (Main script)                                                                                                                         |
|                     | Rawdata/parentalVSReporter/\*.\* (download)                                                                                                                          |
| Fig.1C and S1G      | Fig1C/`fig1CS1G.m` (Main script)                                                                                                                                     |
|                     | Rawdata/western/\*.\* (download)                                                                                                                                     |
| Fig.S1E             | FigS1E/`pulsing_vs_iqr.py` (Main script)                                                                                                                             |
|                     | Rawdata/Workspaces/`130722_SCdyn.csv` (download)                                                                                                                     |
|                     | Rawdata/Workspaces/`130722_Pav.csv` (download)                                                                                                                       |
| Fig.1D and 2A&B     | Fig1D-2AB/`fig1D.m` (Main script)                                                                                                                                    |
|                     | Fig1D-2AB/`fig2A.m` (Main script)                                                                                                                                    |
|                     | Fig1D-2AB/`fig2B.m` (Main script)                                                                                                                                    |
|                     | Rawdata/Workspaces/`site_x.mat` for x = 4, 17, 37, 44, 57, 64(download)                                                                                              |
|                     | Rawdata/Workspaces/`harm_basis_130722_corrected_retracked_all_cleaned_late.mat` (download)                                                                           |
|                     | Rawdata/Workspaces/`fourier_signals_corrected_cleaned_newBTC2.mat` (download)                                                                                        |
| Fig.3 and S2        | Fig3-S2/`fig3A.m` (Main script)                                                                                                                                      |
|                     | Fig3-S2/`fig3B.m` (Main script)                                                                                                                                      |
|                     | Fig3-S2/`fig3C.m` (Main script)                                                                                                                                      |
|                     | Fig3-S2/`fig3D.m` (Main script)                                                                                                                                      |
|                     | Fig3-S2/`figS2A.m` (Main script)                                                                                                                                     |
|                     | Fig3-S2/`figS2B.m` (Main script)                                                                                                                                     |
|                     | Fig3-S2/`figS2C.m` (Main script)                                                                                                                                     |
|                     | Rawdata/Workspaces/`site_x_130722_corrected_retracked_all_cleaned.mat` for x = 17, 57, 64(download)                                                                  |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_newBTC.mat` (download)                                                                                                 |
|                     | Rawdata/Workspaces/`harm_basis_fPCA_5basis_noFGF_newBTC_rot.mat` (download)                                                                                          |
|                     | Rawdata/Workspaces/`harm_basis_50_to_600.mat` (download)                                                                                                             |
|                     | Rawdata/Workspaces/`harm_basis_130722_corrected_retracked_all_cleaned_late.mat` (download)                                                                           |
|                     | Rawdata/Workspaces/`site_4_130722_corrected_retracked_all_paper_cleaned.mat` (download)                                                                              |
|                     | Rawdata/Workspaces/`130722_SCfeat.csv` (download)                                                                                                                    |
| Fig.4, S3 and S5A&B | Fig4-S3-S5AB/`fig4AB.m` (Main script)                                                                                                                                |
|                     | Fig4-S3-S5AB/`fig4CDS3DS5AB.m` (Main script)                                                                                                                         |
|                     | Fig4-S3-S5AB/`figS3A.m` (Main script)                                                                                                                                |
|                     | Fig4-S3-S5AB/`figS3B.m` (Main script)                                                                                                                                |
|                     | Fig4-S3-S5AB/`figS3C_lower.m` (Main script)                                                                                                                          |
|                     | Fig4-S3-S5AB/`figS3C_upper.m` (Main script)                                                                                                                          |
|                     | Rawdata/Workspaces/` site_x_130722_corrected_retracked_all_cleaned.mat` for x = 17, 57, 64(download)                                                                 |
|                     | Rawdata/Workspaces/`harm_basis_fPCA_5basis_noFGF_newBTC_rot.mat` (download)                                                                                          |
|                     | Rawdata/Workspaces/` harm_basis_130722_corrected_retracked_all_cleaned_late_newBTC.mat` (download)                                                                   |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_newBTC.mat` (download)                                                                                                 |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_AKTi.mat` (download)                                                                                                   |
|                     | Rawdata/Workspaces/`scores_puls_corrected_retracked_all_cleaned_newBTC.mat` (download)                                                                               |
|                     | Rawdata/Workspaces/`scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi.mat` (download)                                                                          |
| Fig.5A and S4B&D    | Fig5A-S4BD/`fig5AS4BD.m` (Main script)                                                                                                                               |
|                     | Rawdata/Workspaces/`site_x_04-15-2014_all_paper_cleaned.mat` for x = 1,3,5,7,9,11,14,16,18,20,22,24,25,27,29,31,33,35,38,40,42,44,46,48,49,51,53,55,57,59 (download) |
|                     | Rawdata/Workspaces/`scores_04-15_new.mat` (download)                                                                                                                 |
| Fig.5B and S4C      | Fig5B-S4C/`fig5BS4C.m` (Main script)                                                                                                                                 |
|                     | Rawdata/Workspaces/`site_x_130722_corrected_retracked_all_paper_cleaned.mat` for x = 1,2,4,37,39,40,41,42,44,61,62,64 (download)                                     |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_AKTi.mat` (download)                                                                                                   |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_MEKi.mat` (download)                                                                                                   |
|                     | Rawdata/Workspaces/`scores_early_5basis_noFGF_newBTC.mat` (download)                                                                                                 |
|                     | Rawdata/Workspaces/`scores_puls_corrected_retracked_all_cleaned_newBTC.mat` (download)                                                                               |
|                     | Rawdata/Workspaces/`scores_puls_corrected_retracked_all_cleaned_newBTC_ATKi.mat` (download)                                                                          |
|                     | Rawdata/Workspaces/`scores_puls_corrected_retracked_all_cleaned_newBTC_MEKi.mat` (download)                                                                          |
| Fig.5C              | Fig5C /`fig5C.m` (Main script)                                                                                                                                       |
|                     | Rawdata/Workspaces/`dists_04182014.mat` (download)                                                                                                                   |
| Fig.5D              | Fig5D/`fig5D.m` (Main script)                                                                                                                                        |
|                     | Rawdata/Workspaces/`c_signal_03302014.mat` (download)                                                                                                                |
|                     | Rawdata/Workspaces/`dists_04182014.mat` (download)                                                                                                                   |
|                     | Rawdata/Workspaces/`scores_04182014.mat` (download)                                                                                                                  |
| Fig.6B              | Fig6B/`plot_exampleEKAREVvsF3aN400.m` (Main script)                                                                                                                  |
|                     | Fig6B/`calculate_correlation.m` (Main script)                                                                                                                        |
|                     | Rawdata/dualsensors/\*.\* (download)                                                                                                                                 |
| Fig.6C              | Fig6C/`fig6C.m` (Main script)                                                                                                                                        |
|                     | Rawdata/Workspaces/`140215_SCdyn_rev1.csv` (download)                                                                                                                |
| Fig.6D              | Fig6D/`fig6D.m` (Main script)                                                                                                                                        |
|                     | Rawdata/Workspaces/`scores_03242014.mat` (download)                                                                                                                  |
| Fig.7B              | Fig7/`analysis_median_iqr_rotation.m` (Main script)                                                                                                                  |
|                     | Rawdata/fixedcell/\*.\* (download)                                                                                                                                   |
| Fig. 7C, S7A, S7B   | Fig7/`plot_inhib_effect.m` (Main script)                                                                                                                             |
|                     | Rawdata/fixedcell/\*.\* (download)                                                                                                                                   |
| Fig. 7D, S7E, S7F   | Fig7/`run_sensitivity.m` (Main script)                                                                                                                               |
|                     | Rawdata/fixedcell/\*.\* (download)                                                                                                                                   |