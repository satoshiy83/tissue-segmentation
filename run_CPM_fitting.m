% Function screening parameters for cellular Potts model, returning
% smoothed regions with circularity higher than a given threshold and
% homogeneous as high as possible.
% Written by Satoshi Yamashita.

function result = run_CPM_fitting(partition,data,hint,delegate)
% Function screening parameters for cellular Potts model.
% result = run_CPM_fitting(partition,data,hint,delegate)
% Argument partition is a matrix whose rows and columns represent data
% points and regions at initial state.
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   cpmf_vary_range: (double[1,2]) range to vary temperature and Hamiltonian
%   coefficients.
%   cpmf_vary_zoomF: (double,double[1,2]) zoom factor to narrow the range.
%   cpmf_vary_steps: (double) number of steps to vary parameters.
%   cpmf_vary_convT: (double[1,2]) threshold to check convergence.
%   cpmf_miniR: (double) minimum roundness or regions.
%   cpmf_uppaN: (double) maximum number to update parameters.
%   cpm_iterN: (double) maximum number to try updates.
%   cpm_uplaN: (double) number to update pixel label.
%   cpm_coefficients: (double[1,3]) array of coefficients for area
%   constraint, surface tension, and homogeneity of region.
%   cpm_temperature: (double) tempareture of system.
% Return value is a matrix whose rows and columns represent data points and
% regions smoothed.
% Parameters used to smooth resultant regions were stored in delegate.hint.

%% Initialization.
hint = hint.copy;
range = hint.objectForKey('cpmf_vary_range');
zoomF = hint.objectForKey('cpmf_vary_zoomF');
steps = hint.objectForKey('cpmf_vary_steps');
convT = hint.objectForKey('cpmf_vary_convT');
miniR = hint.objectForKey('cpmf_miniR');
uppaN = hint.objectForKey('cpmf_uppaN');

regiData = SYData(partition);
meter = TSMeter(data,hint);
meter.regiList = regiData;
analyser = TSRegionAnalyser(data);
analyser.regiList = regiData;
analyser.meter = meter;

T = hint.objectForKey('cpm_temperature');
tempArray = lf_temperature(T,range(1),steps,convT(1));
coef = hint.objectForKey('cpm_coefficients');
coefArray = lf_coefficients(coef,range(2),steps,convT(2));

%% Main loop.
while ~(uppaN < 1)
    % try parameters spread in range.
    siz = [tempArray.count,coefArray.count];
    roundness = zeros(siz);
    silhouettes = zeros(siz);
    for i = 1:tempArray.count
        temp = tempArray.objectAtIndex(i);
        hint.setObjectForKey('cpm_temperature',temp.var);
        for j = 1:coefArray.count
            coef = coefArray.objectAtIndex(j);
            hint.setObjectForKey('cpm_coefficients',coef.var);
            qartition = run_CPM_smoothing(partition,data,hint,delegate);
            analyser.regiList = SYData(qartition);
            roun = analyser.roundness_of_regions;
            roundness(i,j) = min(roun);
            silh = analyser.silhouette;
            silhouettes(i,j) = sum(silh(:)) / data.count;
        end
    end
    
    % choose the best parameters in range.
    mask = roundness < miniR;
    if any(~mask(:))
        silhouettes(mask) = nan;
        [~,index] = max(silhouettes(:));
    else
        [~,index] = max(roundness(:));
    end
    
    % update range.
    flag = lf_shall_zoom(siz,index);
    if flag
        range = range .* zoomF;
    end
    [r,c] = ind2sub(siz,index);
    temp = tempArray.objectAtIndex(r);
    tempArray = lf_temperature(temp.var,range(1),steps,convT(1));
    coef = coefArray.objectAtIndex(c);
    coefArray = lf_coefficients(coef.var,range(2),steps,convT(2));
    
    % check convergence.
    if tempArray.count == 1 && coefArray.count == 1
        break;
    end
    
    uppaN = uppaN - 1;
    clc;
end

temp = tempArray.objectAtIndex(1);
hint.setObjectForKey('cpm_temperature',temp.var);
coef = coefArray.objectAtIndex(1);
hint.setObjectForKey('cpm_coefficients',coef.var);
delegate.hint = hint;
qartition = run_CPM_smoothing(partition,data,hint,delegate);

result = qartition;
end

function result = lf_shall_zoom(siz,index)

[r,c] = ind2sub(siz,index);
if siz(1) > 1 && (r == 1 || r == siz(1))
    result = false;
elseif siz(2) > 1 && (c == 1 || c == siz(2))
    result = false;
else
    result = true;
end
end
function result = lf_temperature(T,range,steps,convT)

array = SYArray;

if range == 0
    array.addObject(SYData(T));
    result = array;
    return;
end

while true
    if steps < 2
        array.addObject(SYData(T));
        result = array;
        return;
    elseif range / steps > convT
        break;
    end
    
    steps = steps - 1;
end

M = T + range / 2;
m = T - range / 2;
if m < 0
    m = 0;
end
tArray = m:(M - m) / steps:M;
tArray(tArray <= 0) = [];
for t = tArray
    array.addObject(SYData(t));
end
result = array;
end
function result = lf_coefficients(coef,range,steps,convT)

array = SYArray;

if range == 0
    array.addObject(SYData(coef));
    result = array;
    return;
end

while true
    if steps < 2
        array.addObject(SYData(coef));
        result = array;
        return;
    elseif range / steps > convT
        break;
    end
    
    steps = steps - 1;
end

M = range / 2;
m = -range / 2;
area_constraint = coef(1);
surface_tension = coef(2);
silhouette = coef(3);
for p = m:range / steps:M
    k = 2 ^ p;
    array.addObject(SYData([area_constraint, ...
        surface_tension * k, ...
        silhouette / k]));
end
result = array;
end
