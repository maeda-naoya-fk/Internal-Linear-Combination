# ILC(Internal Linear Combination)

***Nonparametric component separation for Cosmic Microwave Background observations in harmonic space***


## Description

ILC(Internal Linear Combination) is the foreground removal method.  
You can subtract foreground and noise from observation imaps using this method, and get clean CMB map.  
You should note that this method is perforemed in **harmonic space**.  
For detail, see paper 'A high resolution foreground cleaned CMB map from WMAP'(T. Max, O. Angelica et al., astro-ph/0302496)

## Usage

You can just run this code, and get result.
If you want to change parameter, Please note that maximum of l is 750 because lmax of beam_transfer_function wmap team produced is 750.

## Requirement

- numpy(=1.18.1)
- healpy(=1.12.10)
- matplotlib(=3.2.1)

## Contributors

- Eiichiro Komatsu (MPA, KavliIPMU)


