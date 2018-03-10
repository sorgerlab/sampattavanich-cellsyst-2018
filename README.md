# Sampattavanich et al., Cell Systems (2018)

**Source code and example plots for the article "Encoding growth factor
identity in the temporal dynamics of a transcription factor under
combinatorial regulation" by Sampattavanich et al.**

Each piece of source code is provided in a folder containing matlab
scripts and all related functions. The main function to generate related
figures are highlighted in **bold** in the following table. Readme file
can also be seen in the accompanied PDF file, also showing example
plots. CellProfiler project files (`.cpproj`) are also provided for users
who are interested to see our pipelines for image segmentation.

To run this code, users must download related source files and put these
in the \\rawdata folder placed at the top-most level of this git
repository folder.

| Figures           | Related files                                         |
|-------------------|-------------------------------------------------------|
| Fig.1B            | `createParentalVSReporter.m` (Main code)              |
|                   | `genFoxO3.m`                                          |
|                   | `plotParentalVSReporter.m`                            |
|                   | `errorbarxy.m` by Qi An (matlabcentral)               |
|                   | `createParentalVSReporter-readme.pdf`                 |
|                   | Rawdata/`combineddata06012014.mat` (download)         |
|                   | Rawdata/`analysisPipe06012014-parental-withBG.cpproj` |
|                   | Rawdata/`analysisPipe06012014-reporter-withBG.cpproj` |
| Fig.S1B and S1C   | `comparedataat15min.m` (Main code)                    |
|                   | `Comparedataat15min-readme.pdf`                       |
|                   | Rawdata/`combineddata06012014.mat` (download)         |
|                   | Rawdata/`analysisPipe06012014-parental-withBG.cpproj` |
|                   | Rawdata/`analysisPipe06012014-reporter-withBG.cpproj` |
| Fig.7B            | ‘analysis_median_iqr_rotation.m’ (Main code)          |
|                   | Rawdata/fixedcell/\*.\* (download)                    |
| Fig. 7C, S7A, S7B | `plot_inhib_effect.m` (Main code)                     |
|                   | Rawdata/fixedcell/\*.\* (download)                    |
| Fig. 7D, S7E, S7F | `run_sensitivity.m` (Main code)                       |
|                   | Rawdata/fixedcell/\*.\* (download)                    |