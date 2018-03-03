# Use Case 02

Pipeline example to analyze fNIRS data exported in *.nirs format without trigger markers information

### Input:
  - Insert markers in *.nirs files according to the experimental protocol (A/B)
  - Load all *.nirs files organized in the group-subject folder hierarchy
  - Create demographics table with information about groups and subjects
  
### Preprocessing:
  - Set duration of the task to 20s for all subjects and repetitions
  - Remove auxiliary (extra) trigger marker from the loaded data
  - Truncate time series (30s prior first onset, 30s after last onset + duration)
  - Compute the hemodynamic states. Note: unfiltered for AR-IRLS modelling
  
### Statistical Analysis:
  - Create the GLM for each subject based on AR-IRLS and the canonical hrf
  - Generate individual results for p < 0.05 as corrected for Bonferroni
  - Create a mixed effects model for the group level analysis
  - Generate the group results for p < 0.05 as corrected for Bonferroni
  
### Data Visualization:
  - Band-pass filter the hemodynamic data for visualization purposes
  - Generate a block average for each subject to inspect activation

.
**Important note**: the script is primarily built upon the nirs-toolbox by Ted Huppert et al. [1]
[1]: https://bitbucket.org/huppertt/nirs-toolbox/wiki/Home

**Bugs reports and suggestions for improvements would be very appreciated.**