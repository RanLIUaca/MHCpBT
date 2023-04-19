# MHCpBT
A Bayesian Approach to Estimate MHC-Peptide Binding Threshold

## Installation
You can install the released version of MHCpBT from GitHub with:
```
devtools::install_github("RanLIUaca/MHCpBT")
library(MHCpBT)
```

## Prior settings
We set the non-informative prior distributions as the default settings for users: $\boldsymbol{\gamma}$ and $\boldsymbol{\eta_j}$ are two vector whose entries are all one. $\boldsymbol{\mu_0}$ is a zero vector, and $\boldsymbol{\Sigma_0}$ is a diagonal matrix, $10\boldsymbol{I}$. $\alpha_0$ and $T_0$ are both set to one. 

## Input
data1: an array whose elements are strings representing the peptide sequences. Each element of the array must have the same length.
Y: a numeric vector representing the binding scores of the peptides in data1. The length of Y must be equal to the length of data1.
motif_len: an integer representing the length of the binding motif. This is a required input.
dict: a named character vector containing the amino acid letters and their corresponding codes. The default value is c(A=1,C=2,D=3,E=4,F=5,G=6,H=7,I=8,K=9,L=10,M=11,N=12,P=13,Q=14,R=15,S=16,T=17,V=18,W=19,Y=20).
burn_in: an integer representing the number of iterations to discard as burn-in. The default value is 2000.
end_times: an integer representing the total number of iterations to run the Gibbs sampler. The default value is 3000.
result_path: a character string representing the directory where the results will be stored. The default value is "./".
only_get_func: a logical value that specifies if only the functions for getting the results should be returned. The default value is FALSE.

## Usage
The main function in MHCpBT is threshold_gibbs, which estimates the MHC-peptide binding threshold using Gibbs sampling.
```
threshold_gibbs(data1, Y, motif_len, dict, burn_in=2000, end_times=3000, result_path,
                          only_get_func = F,input_prior = F,prior = 0)
```

## Example
This is a toy example for estimating the binding threshold for a specific allele:
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
-   Ran, L. et al. (2023+), “A Bayesian Approach to Estimate MHC-Peptide Binding Threshold,” Working Paper.

