# tissue-segmentation
This project provides functions for automatic segmentation of epithelial tissue during morphogenesis.

## Installation
Download files and put them in a folder with a suitable name. Go to Matlab command line, enter "addpath" + a full path of the folder, and enter "savepath".

## Requirement
This project requires no Matlab toolbox, but a custom framework of objective classes SYObject family which is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3625631.svg)](https://doi.org/10.5281/zenodo.3625631).

## Usage
Assume that *Q* be a 3D matrix of 100 times 100 times 200 rows and columns and slices. To segment the 100 times 100 rows and columns based on L^1 norm between the 200D vectors, prepare a dictionary holding parameters and run region growing, label propagation, and cellular Potts model algorithms.

The parameters are given in a dictionary.
```
hint = SYDictionary;
hint.setObjectForKey('property_depth',1);
hint.setObjectForKey('time_length',200);
hint.setObjectForKey('neigN',8);
hint.setObjectForKey('metric_base','time_integral_absolute');
hint.setObjectForKey('metric_exte','mean');
hint.setObjectForKey('collector','centroids');
hint.setObjectForKey('rg_regiN',6);
hint.setObjectForKey('rg_inseN',50);
hint.setObjectForKey('rg_upseN',20000);
hint.setObjectForKey('rg_cyclN',20);
hint.setObjectForKey('lp_regiN',6);
hint.setObjectForKey('lp_upthN',20000);
hint.setObjectForKey('lp_upcmN',20000);
hint.setObjectForKey('lp_iterN',1);
hint.setObjectForKey('lp_uplaN',20000);
hint.setObjectForKey('cpmf_vary_range',[2,4]);
hint.setObjectForKey('cpmf_vary_zoomF',[0.5,0.5]);
hint.setObjectForKey('cpmf_vary_steps',4);
hint.setObjectForKey('cpmf_vary_convT',[0.1,0.1]);
hint.setObjectForKey('cpmf_miniR',0.45);
hint.setObjectForKey('cpmf_uppaN',20000);
hint.setObjectForKey('cpm_iterN',20000);
hint.setObjectForKey('cpm_uplaN',50);
hint.setObjectForKey('cpm_coefficients',[0.001,10,1]);
hint.setObjectForKey('cpm_temperature',1);
```
For each paramter, look up comments inside functions.
