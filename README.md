# Mada-Ct-Distribute

This repository contains all scripts and data to support a recent submission (Andriamandimby and Brook et al. 2021) to Epidemics' special issue on COVID-19 dynamics in lower- and middle- income countries. Each subfolder of the repository contains scripts for largely self-contained analyses. These are organized as follows:

* *epinow2-ipm*: This folder contains COVID-19 case count data from March - September 2020 at the National level and across two main administrative regions, derived from RT-qPCR testing in the Virology Unit laboratory at Institut Pasteur of Madagascar. The linelist of cases ('ipm-case-dat.csv') is then called by three R-scripts "estim-*" which use the opensource EpiNow2 package to estimate Rt and the epidemic growth rate from each region. The output .Rdata files are also included in the folder. A summary script "extract-save-growth-rt.R" compiles these Rdata files into a final csv which is referenced in subsequent analyses. The .sh files used to run this script on the UC Berkeley savio cluster are also included.

* *epinow-reported*: This folder contains Rt estimates from national data for Madagascar ("rt_ests.csv") computed using the pipeline available at: https://github.com/labmetcalf/mada-rt. Also included is an R script that converts these Rt estimates to COVID-19 growth rates and compiles them with estimates derived from the IPM data referenced above.

* *TC-dat-up*: This folder contains raw data and processing script from tissue culture inoculations used to correct all Ct values to TaqPath N gene values.

* *common_data*: This folder contains all major Ct datasets used in genrealized additive modeling, as well as the Ct correction parameters derived from TC analysis.

* *Ct-longitudinal-gams*: Contains data, scripts, and compilation script for longitudinal GAMs of population-level Ct by region (Fig. 2A, main text, Fig. S2, Table S4). Also includes .sh scripts to run these analyses on Savio or any other remote cluster. Note that actual fitted GAMs are not included due to large filesize but should be reproducible from data and scripts included (the National dataset will take a long time, ~24hours, to run).

* *individual-Ct*: Contains individual-Ct trajectory data and GAM scripts, as well as viral kinetics model fitting script adopted from Hay et al. 2020 (see [here](https://github.com/jameshay218/virosolver_paper) for original script). Note that actual fitted GAMs are not included due to large filesize but should be reproducible from data and scripts included. Output of these analyses is summarized in Fig 2B-D main text and Table S5. Folder also includes .sh scripts to run these analyses on Savio or any other remote cluster.

* *common_pars*: Contains fitted viral kinetics parameters derived from mechanistic viral kinetics model fits to the individual trajectory GAM output.

* *Ct-sym-asym*: Contains data and script for cross-sectional population-level GAM querying the effect of symptom status (independent of duration of infection) on Ct. Note that actual fitted GAM is not included due to large filesize but should be reproducible from data and scripts included. Output is summarized in Table S6 of the manuscript.

* *viro-gp-clust*: Contains scripts, data, and .sh savio scripts to run Gaussian process model fit to cross-sectional Ct distributions across all three regions. Script is adopted from [Hay et al. 2020](https://github.com/jameshay218/virosolver_paper).

* *viro-seir-clust*: Contains scripts, data, and .sh savio scripts to run SEIR model fit to cross-sectional Ct distributions across all three regions. Script is adopted from [Hay et al. 2020](https://github.com/jameshay218/virosolver_paper).

* *fig-plots*: Contains scripts and resulting figures to produce all main text and supplementary figures in the manuscript. Note that GAMs will need to be run in advance to produce fitted GAMs in scripts where these are called as .Rdata files. GIS shapefiles for Madagascar (Fig. 1) are not included here due to large filesize but are available for public download as MDG_adm0 (boundary) and MDG_adm3 (third-level administrative boundaries) [here](https://earthworks.stanford.edu/catalog.html?f%5Bdc_format_s%5D%5B%5D=Shapefile&f%5Bdc_subject_sm%5D%5B%5D=Boundaries&f%5Bdct_spatial_sm%5D%5B%5D=Madagascar&per_page=20&sort=dc_publisher_sort+asc%2C+dc_title_sort+asc).


