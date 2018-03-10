**Sampattavanich et al., Cell Systems (2018)**

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

<table>
<thead>
<tr class="header">
<th>Figures</th>
<th>Related files</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Fig.1B</td>
<td><p><strong>`createParentalVSReporter.m`</strong></p>
<p>`genFoxO3.m`</p>
<p>`plotParentalVSReporter.m`</p>
<p>`errorbarxy.m` by Qi An (matlabcentral)</p>
<p>`createParentalVSReporter-readme.pdf`</p>
<p>Rawdata/`combineddata06012014.mat` (<strong>download</strong>)</p>
<p>Rawdata/`analysisPipe06012014-parental-withBG.cpproj`</p>
<p>Rawdata/`analysisPipe06012014-reporter-withBG.cpproj`</p></td>
</tr>
<tr class="even">
<td>Fig.S1B and S1C</td>
<td><p><strong>`comparedataat15min.m`</strong></p>
<p>`Comparedataat15min-readme.pdf`</p>
<p>Rawdata/`combineddata06012014.mat` (<strong>download</strong>)</p>
<p>Rawdata/`analysisPipe06012014-parental-withBG.cpproj`</p>
<p>Rawdata/`analysisPipe06012014-reporter-withBG.cpproj`</p></td>
</tr>
<tr class="odd">
<td>Figure 7C, S7A, S7B</td>
<td><p><strong>`plot_inhib_effect.m`</strong></p>
<p>Rawdata/fixedcell/*.* (<strong>download</strong>)</p></td>
</tr>
<tr class="even">
<td>Figure 7D, S7E, S7F</td>
<td><p><strong>`run_sensitivity.m`</strong></p>
<p>Rawdata/fixedcell/*.* (<strong>download</strong>)</p></td>
</tr>
</tbody>
</table>
