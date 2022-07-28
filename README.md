**A Kernel dIfferentiability correlation-based NOn-negative Matrix factorization algorithm using Kullback-Leibler divergence loss Optimization (KINOMO)**

We propose Kernel dIfferentiability correlation-based NOn-negative Matrix factorization algorithm using Kullback-Leibler divergence loss Optimization (KINOMO), a semi-supervised NMF model that is robust to noises and also uses ‘prior’ biological knowledge for better refinement.

**Instructions for running KINOMO**
Please note that KINOMO needs to be implemented sample-wise and not on the integrated object (if at all it exists).

**Accessing the KINOMO repository**

1. The KINOMO repsository can be accessed via https://github.com/IzarLab/KINOMO.git
2. Download the KINOMO repository into your local machine.


**Running KINOMO**
1. Detailed step-by-step instructions/help (in R) for running KINOMO is specified in the **'kinomo_run.R'**. Please use this script.
2. If there is any error with respect to accessing specific sub-scripts/modules, please provide necessary path-to-file information.

**Packages required**
dplyr
Seurat
ggplot2
gplots
purrr
cowplot
stringr
NMF
pkgmaker
cluster

**Input file**
Input file is the gene expression data (raw counts) converted as a Seurat object.

**Estimating the factorization rank**
1. A critical parameter in KINOMO is the factorization rank r. It defines the number of metagenes used to approximate the target matrix. 
2. Given a NMF method and the target matrix, a common way of deciding on r is to try different values, compute some quality measure of the results, and choose the best value according to this quality criteria.
3. Several approaches have then been proposed to choose the optimal value of r. For example, (Brunet2004) proposed to take the first value of r for which the cophenetic coefficient #starts decreasing, (Hutchins2008) suggested to choose the first value where the RSS curve presents an inflection point, and (Frigyesi2008) considered the smallest value at which #the decrease in the RSS is lower than the decrease of the RSS obtained from random data.

**Performing factorization using estimated rank**
1. Based on the consensus rank estimation, factorization is done.

**Metagene identification**
1. top50/100/200/300 meta-genes are estimated.
2. Meta-gene signatures are plotted
3. Meta-gene heatmaps are plotted

**After the script successfully runs, the following files would be created in the current directory:**
1. Rank plot
2. Meta-gene signature UMAPS
3. Meta-gene heatmaps
4. Associated .csv files
5. Associated .rds files
