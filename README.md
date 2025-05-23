# MHCpBT
MHCpBT is an R package for estimating the MHC-peptide binding threshold using Gibbs sampling. This package is designed to identify an appropriate binding score that differentiates strong and weak binding peptides. The package takes as input peptide sequences and corresponding binding scores and outputs the estimated binding threshold. This package is also suitable for other biological binding sequences.

## Installation
You can install the released version of MHCpBT from GitHub with:
```
devtools::install_github("RanLIUaca/MHCpBT")
library(MHCpBT)
```

## Prior settings
We set the non-informative prior distributions as the default settings for users: $\boldsymbol{\gamma}$ and $\boldsymbol{\eta_j}$ are two vector whose entries are all one. $\boldsymbol{\mu_0}$ is a zero vector, and $\boldsymbol{\Sigma_0}$ is a diagonal matrix, $10\boldsymbol{I}$. $\alpha_0$ and $T_0$ are both set to one. 

## Input
* data1: an array whose elements are strings representing the peptide sequences.
* Y: a numeric vector representing the binding scores of the peptides in data1. The length of Y must be equal to the length of data1.
* motif_len: an integer representing the length of the binding motif. This is a required input.  Usually, it is 9.
* dict: a named character vector containing the amino acid letters and their corresponding codes. 
* burn_in: an integer representing the number of iterations to discard as burn-in. The default value is 2000.
* end_times: an integer representing the total number of iterations to run the Gibbs sampler. The default value is 3000.
* result_path: a character string representing the directory where the results will be stored. The default value is "./".
* only_get_func: a logical value that specifies if only the functions for getting the results should be returned. The default value is FALSE.

## Output
The output of this code includes:

- Several diagnostic plots:
  - A plot of the final peptide logo
  - A plot of the joint posterior trace
  - A plot of the coefficient trace
  - A plot of the MSE trace
  - A plot of the real vs estimated scores
  - A plot of the correlation between the true and estimated scores
  - A plot of the estimated binding probabilities
- CSV files containing results of the analysis:
  - `temp_X.csv`: temporary X matrix
  - `temp_beta.csv`: temporary beta coefficients
  - `result_threshold.csv`: estimated binding threshold
  - `d_trace.csv`: trace of the estimated binding threshold
  - `log_joint_post_path.csv`: joint posterior distribution trace
  - `coe_path.csv`: coefficient trace
  - `MSE_path.csv`: MSE trace
  - `result_X.csv`: $\boldsymbol{X}$ matrix 
  - `result_A.csv`: estimated $\boldsymbol{A}$ matrix
  - `result_beta.csv`: estimated $\boldsymbol{\beta}$ 
  - `result_theta.csv`: estimated $\boldsymbol{\Theta}$ 
  - `result_sigma2.csv`: estimated $\boldsymbol{\sigma^2}$
  - `result_theta_0.csv`: estimated $\boldsymbol{\theta_0}$
  - `result_score.csv`: estimated binding scores
  - `final_result.csv`: final data frame including peptide sequences, binding scores, estimated binding probabilities, and estimated binding labels
  - `result_coefficient.csv`: correlation coefficients between the true and estimated binding scores

## Usage
The main function in MHCpBT is threshold_gibbs, which estimates the MHC-peptide binding threshold using Gibbs sampling.
```
threshold_gibbs(data1, Y, motif_len, dict, burn_in=2000, end_times=3000, result_path,
                          only_get_func = F,input_prior = F,prior = 0)
```

## Example
This is a toy example for estimating the binding threshold for a specific allele. We shall use the data contained in the package.
```
library('MHCpBT')
library('stringr')
set.seed(2020)

# settings
motif_len = 9
allele_name = 'H-2-Ld'
burn_in = 1
end_times = 10
result_path = 'result/'
if(dir.exists(result_path)==0){dir.create(result_path)}

# read data
identity_threshold = 80
data_path = system.file("extdata", paste0('ba_', allele_name, identity_threshold,'.csv'), package = "MHCpBT")
data1 = read.csv(data_path,stringsAsFactors = F)
data1 = data1[,1:2]
colnames(data1) =  c('pep','score')
strl = str_length(data1$pep)
data1 = data1[which(strl>=motif_len),]
Y = data1$score
data1 = data1$pep

# AA
dict = 'ACDEFGHIKLMNPQRSTVWY'
# DNA
# dict = 'ATGC'
tmp_index = 1:str_length(dict)
dict = str_sub(dict,tmp_index,tmp_index)
digit_dict = as.character(1:length(dict))
names(digit_dict) = dict
names(dict) = digit_dict 

result = threshold_gibbs(data1 = data1, Y = Y, motif_len, dict = dict,
                         burn_in=burn_in, end_times=end_times, 
                         result_path = result_path)
```

## Reference
-   Liu, R., Hu, Y., Huang, J., & Fan, X. (2023). A Bayesian approach to estimate MHC-peptide binding threshold. Briefings in Bioinformatics, 24(4), bbad208. https://doi.org/10.1093/bib/bbad208

