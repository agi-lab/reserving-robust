# reserving-robust
R code for academic paper written by Avanzi, Lavender, Taylor and Wong: *Detection and treatment of outliers in multivariate robust loss reserving*.

This code has been arranged so that it can be run sequentially, where appropriate outliers and results of the detection methods can be viewed as each line is run. For both scripts, the code applies the *robust bivariate chain-ladder* technique (Verdonek and Van Wouwe, 2011) to the different outlier detection and adjustment methods explained in the paper.

## Bivariate Code and dataset
This file contains the code for generating Sections 2 and 3 of the academic paper. The data used in this section *(ShBaMe12 Robust Bivariate CL Data.csv)* is taken from Shi, Basu, and Meyers (2012), containing the incremental claims for Personal and Commerical auto insurance from a major U.S. property-casualty insurer from 1988-1997.

## Trivariate Code and dataset
This file contains the code for generating Sections 4 of the academic paper, showing the outlier detection and treatment for trivariate data. The data used in this section *(Trivariate Data.csv)* is synthetic and aims to simulate incremental claims from the same three lines of businesses analysed in the paper.

