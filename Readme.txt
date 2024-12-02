The R codes and the data files used in Srivastava, Sengupta, and Ghosh (2024), "Quantifying the rhythm of ice accumulation in East Antarctica with varying temperature" is given in this repository. The supplementary material of this article is also available here for reference.

The three data sets are publicly available from the sites mentioned in the references. 
The AICC2012 age scales for Lake Vostok and EPICA Dome C (aicc2012icecore-data.xls) were downloaded from https://www.ncei.noaa.gov/access/paleo-search/study/15076. 

Data files

Dome_Fuji_downloaded.csv: This file contains the downloaded data.
Vostok_downloaded.csv: This file contains the downloaded data. The age in this file is AICC2012 ages.
EPICA_downloaded.csv: This file contains the downloaded data. The age in this file is AICC2012 ages.     

Dome_Fuji_data.csv, Vostok_data.csv, EPICA_data.csv : This file contains log(AAR) (y), temperature (x) and age(z). These columns are derived from the downloaded data.


Figure2_Fuji.R, Figure2_Vostok.R and Figure2_EPICA.R files are R script files that read the data files and produce plots of Figure 2.

The R script file Estimation_functions.R consists of several R functions. 
The function LSE computes the least square estimates of parameters. 
The function Estimates computes the smooth least square estimates of g and other parameters.
The function Estimates_Bootstrap generate the model based bootstrap samples and compute the uncertaininity quantities related to the parameters.

The R script files Estimation_plots_Dome_Fuji.R, Estimation_plots_Vostok.R, Estimation_plots_EPICA.R uses the data (log(AAR) (y), temperature (x) and age (z)) and Estimation_functions.R for computation of parameters, uncertainity estimates and plots of Figure 3 and 4.
