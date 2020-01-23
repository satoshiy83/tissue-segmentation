# tissue-segmentation
This project provides functions for automatic segmentation of epithelial tissue during morphogenesis.

## Installation
Download files and put them in a folder with a suitable name. Go to Matlab command line, enter "addpath" + a full path of the folder, and enter "savepath".

## Requirement
This project requires no Matlab toolbox, but a custom framework of objective classes SYObject family which is available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3625631.svg)](https://doi.org/10.5281/zenodo.3625631).

## Example
Assume that *Q* be a 3D matrix of 100 times 100 times 200 rows and columns and slices. To segment the 100 times 100 rows and columns into 6 regions based on L^1 norm between the 200D vectors, prepare a dictionary holding parameters and run region growing, label propagation, and cellular Potts model algorithms.

The parameters are given in a dictionary.
```
hint = SYDictionary;
hint.setObjectForKey('property_depth',200);
hint.setObjectForKey('time_length',1);
hint.setObjectForKey('neigN',8);
hint.setObjectForKey('metric_base','absolute');
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

Convert the 3D matrix to a model object.
```
dataMap = TSDataMap(Q,hint);
```

Divide dataMap 50 times by *run_region_growing*(), and save the result into file named *result.rg.mat*.
```
delegate = TSRegionGrowingDelegate;
stack = run_region_growing(dataMap,hint,delegate);
data = SYData(stack);
data.writeToFile('result.rg.mat');
```
You can later read the data in a file written by SYData instanc with an SYData method *initWithContentsOfFile*().
```
data = SYData;
data.initWithContentsOfFile('result.rg.mat');
stack = data.var;
```

Integrate the 50 trials of region growing　by *run_cm_thresholding_lp*(), and save the result into *result.lp.mat* and *result.lp.png*.
```
partition = run_cm_thresholding_lp(stack,hint);
data = SYData(partition);
data.writeToFile(result.lp.mat');
image = regiList2image(partition,dataMap,[100,100]);
image = ss_convert_stack_to_hue(image);
image.writeToFile('result.lp.png',false);
```

To smooth boundaries, first screen paramteres for the cellular Potts model by *run_CPM_fitting*(), and save the parameters to *hint.mat* and *hint.txt*.
```
delegate = CPMDelegate;
run_CPM_fitting(partition,dataMap,hint,delegate);
hint = delegate.hint;
data = hint.data;
data.writeToFile('hint.mat');
txt = hint.description;
fid = fopen('hint.txt','w');
fprintf(fid,txt);
fclose(fid);
```
You can later retrieve the parameter dictionary from the data.
```
data = SYData;
data.initWithContentsOfFile('hint.mat');
hint = SYDictionary;
hint.initWithData(data);
```
With the screened parameters, iterate the boundary smoothing by *run_CPM_smoothing*() 50 times, integrate its results by the label propagation, and save its result to *result.cpm.mat* and *result.cpm.png*.
```
k = hint.objectForKey('lp_regiN');
n = hint.objectForKey('rg_inseN');
stack = zeros(dataMap.count,k,n);
for i = 1:n
  delegate = CPMDelegate;
  qartition = run_CPM_smoothing(partition,dataMap,hint,delegate);
  stack(:,1:size(qartition,2),i) = qartition;
end
partition = run_cm_thresholding_lp(stack,hint);
data = SYData(partition);
data.writeToFile('result.cpm.mat');
image = regiList2image(partition,dataMap,[100,100]);
image = ss_convert_stack_to_hue(image);
image.writeToFile('result.cpm.png',false);
```
Note that the partition is used recursively inside the for loop and so the result of run_CPM_smoothing() is passed to **q**artition.

To calculate silhouette values, use a supporting object TSRegionAnalyser.
```
regiList = SYData(partition);
meter = TSMeter(dataMap,hint);
analyser = TSRegionAnalyser(dataMap);
analyser.meter = meter;
meter.regiList = regiList;
analyser.regiList = regiList;
silhouette = analyser.silhouette;
data = SYData(silhouette);
data.writeToFile('silhouette.cpm.mat');
image = ss_draw_silhouette_map(partition,dataMap,hint);
image = ss_color_blue_red(image);
image.writeToFile('silhouette.png',false);
```

To make a control, use *run_nonsense_segmentation*().
```
hint.setObjectForKey('ns_regiN',6);
hint.setObjectForKey('ns_bootN',20000);
hint.setObjectForKey('ns_insoC',1);
hint.setObjectForKey('ns_fh_score','lf_silhouette');
hint.setObjectForKey('hg_stepN',100);
hint.setObjectForKey('hg_sortD','descend');
scores = run_nonsense_segmentation(dataMap,hint);
data = SYData(scores);
data.writeToFile('silhouette.ctrl.mat');
hg = histogram_scores(scores,hint);
data = SYData(hg);
data.writeToFile('silhouette_histogram.ctrl.mat');
```
