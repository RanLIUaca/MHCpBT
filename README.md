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

