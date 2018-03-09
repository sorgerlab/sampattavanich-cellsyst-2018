**Comparedataat15min.m, relates to supplementary Fig.S1B and S1C**

Required input file: combineddata06012014.mat

Cell lines: parental 184A1 and 184A1 with F3aN400 sensor

Treatment Condition: Cells were serum starved for 6 hours, then treated
with combination of drugs and ligands for 0, 15, 45, 90, 135 minutes.
Parental 184A1 cells were immunostained with FoxO3 antibody. Reporter
184A1 cells with F3aN400 sensor were fixed and imaged.

Cell Profiler (Ver:2.1.1) code for image segmentation and data
extraction

-   Rawdata/analysisPipe06012014-parental-withBG.cpproj
-   Rawdata/analysisPipe06012014-reporter-withBG.cpproj

Contains: 2 struct variables: parental and reporter

+-----------------------------------+-----------------------------------+
| reporter =                        | parental =                        |
|                                   |                                   |
| struct with fields:               | struct with fields:               |
|                                   |                                   |
| Time0: {8×12 cell}                | Time0: {8×12 cell}                |
|                                   |                                   |
| Time15: {8×12 cell}               | Time15: {8×12 cell}               |
|                                   |                                   |
| Time45: {8×12 cell}               | Time45: {8×12 cell}               |
|                                   |                                   |
| Time90: {8×12 cell}               | Time90: {8×12 cell}               |
|                                   |                                   |
| Time120: {8×12 cell}              | Time120: {8×12 cell}              |
+===================================+===================================+
| Dimensions of each field follow the platemap above.                   |
+-----------------------------------+-----------------------------------+
|                            Common data fields                         |
+-----------------------------------+-----------------------------------+
| \'nuc\_area\'                     | \'extendedcyto\_area\'            |
|                                   |                                   |
| \'nuc\_coordX\'                   | \'extendedcyto\_integratedDAPI\'  |
|                                   |                                   |
| \'nuc\_coordY\'                   | \'extendedcyto\_integratedFoxO3\' |
|                                   |                                   |
| \'nuc\_formfactor\'               | \'extendedcyto\_integratedECDGree |
|                                   | n\'                               |
| \'nuc\_integratedDAPI\'           |                                   |
|                                   | \'extendedcyto\_meanDAPI\'        |
| \'nuc\_integratedFoxO3\'          |                                   |
|                                   | \'extendedcyto\_meanFoxO3\'       |
| \'nuc\_integratedWCDGreen\'       |                                   |
|                                   | \'extendedcyto\_meanWCDGreen\'    |
| \'nuc\_meanDAPI\'                 |                                   |
|                                   | \'smallcyto\_area\'               |
| \'nuc\_meanFoxO3\'                |                                   |
|                                   | \'smallcyto\_integratedDAPI\'     |
| \'nuc\_meanWCDGreen\'             |                                   |
|                                   | \'smallcyto\_integratedFoxO3\'    |
| \'nuc\_normmeanFoxO3\'            |                                   |
|                                   | \'smallcyto\_integratedWCDGreen\' |
| \'log10CoverN\_extended\_norm\'   |                                   |
|                                   | \'smallcyto\_meanDAPI\'           |
| \'log10CoverN\_extended\'         |                                   |
|                                   | \'smallcyto\_meanFoxO3\'          |
| \'log10CoverN\_4pixel\_norm\'     |                                   |
|                                   | \'smallcyto\_meanWCDGreen\'       |
| \'log10CoverN\_4pixel\'           |                                   |
+-----------------------------------+-----------------------------------+

