% TSRegionAnalyser: objective class to analyse regions.
% Written by Satoshi Yamahsita.

classdef TSRegionAnalyser < SYObject
properties
    data = nan; % (TSDataMap) given at initialization.
%     hint = nan; % (SYDictionary) given at initialization.
    
    property_depth = []; % given in the hint.
    time_length = []; % given in the hint.
    regiList = []; % (SYData) given after init.
    meter = nan; % (TSMeter) given after init.
    regiN = []; % given in the hint.
end

methods
function obj = TSRegionAnalyser(dataMap)
% Tissue segmentation class to analyse partition.
% obj = TSRegionAnalyser(data,hint)
% Argument data is a TSDataMap instance.
    if nargin < 1
        return;
    end
    
    obj.initWithData(dataMap);
end
function obj = initWithData(obj,dataMap)
% Initialization method with TSDataMap instance.
% obj = initWithData(obj,data)
% Argument data is a TSDataMap instance.
    obj.data = dataMap;
%     obj.hint = hint;
    
    obj.property_depth = dataMap.property_depth;
    obj.time_length = dataMap.time_length;
%     obj.regiN = hint.objectForKey('regiN');
end

function dest = copy(obj,dest)
% Method to make a copy.
    if nargin < 2
        dest = TSRegionAnalyser;
    end
    copy@SYObject(obj,dest);
    
    dest.data = obj.data;
%     dest.hint = obj.hint;
    dest.regiList = obj.regiList;
    dest.meter = obj.meter;
    dest.regiN = obj.regiN;
end

function set.regiList(obj,newList)
    obj.regiList = newList;
    obj.setRegiList(newList);
end
function setRegiList(obj,newList)
    obj.regiN = size(newList.var,2);
end

function result = size_regions(obj)
% Method returning sizes of each region.
% result = size_regions(obj)
% Return value is a row vector.
    result = sum(obj.regiList.var,1);
end
function result = average_regions(obj)
% Method returning average value of each region.
% result = average_regions(obj)
% Return value is a matrix in which row vectors represent the averages.
    aver = zeros(obj.regiN,obj.property_depth * obj.time_length);
    
    for i = 1:obj.regiN
        regi = obj.data.dataList(obj.regiList.var(:,i),3:end);
        aver(i,:) = mean(regi,1);
    end
    
    result = aver;
end
function result = standard_deviation_of_regions(obj)
% Method returning standard deviation of each region.
% result = standard_deviation_of_regions(obj)
% Return value is a column vector.
% Class instance must have an instance variable (TSMeter) meter.
    dev = zeros(obj.regiN,1);
    
    aver = obj.average_regions;
    sizes = obj.size_regions;
    dict = obj.meter.hint.copy;
    obj.meter.hint.setObjectForKey('TSMeterPointConversionSkip',true);
    
    for i = 1:obj.regiN
        array = obj.meter.distanceArray(aver(i,:));
        array = array(obj.regiList.var(:,i));
        dev(i) = sqrt(sum(array .^ 2) / (sizes(i) - 1));
    end
    
    obj.meter.hint = dict;
    result = dev;
end

function result = rim_points(obj)
% Method returning points on regions rim.
% result = rim_points(obj)
% Return value is a logical column vector positive for rim points.
    array = zeros(obj.data.count,1,'logical');
    for i = 1:obj.data.count
        neighbors = obj.data.neighborsOfPoint(i);
        neighbors = obj.regiList.var(neighbors,:);
        if sum(obj.regiList.var(i,:) | any(neighbors,1)) > 1
            array(i) = 1;
        end
    end
    
    result = array;
end
function result = perimeter_length(obj)
% Method returning lengths of region perimeters.
% result = perimeter_length(obj)
% Return value is a row vector containing lengths of region perimeters.
    array = zeros(1,obj.regiN);
    siz = obj.data.frameSize;
    bitmap = zeros(siz,'logical');
    citmap = zeros(siz + 2,'logical');
    indices = 2:(siz(1) + 1);
    jndices = 2:(siz(2) + 1);
    m = [-1,0; 0,-1; 1,0; 0,1];
    for i = 1:obj.regiN
        bitmap(:) = false;
        kndices = obj.regiList.var(:,i);
        r = obj.data.dataList(kndices,2);
        c = obj.data.dataList(kndices,1);
        kndices = sub2ind(siz,r,c);
        bitmap(kndices) = true;
        
        citmap(:) = false;
        for j = 1:size(m,1)
            citmap(indices + m(j,1),jndices + m(j,2)) = ...
                citmap(indices + m(j,1),jndices + m(j,2)) | bitmap;
        end
        citmap(indices,jndices) = citmap(indices,jndices) & ~bitmap;
        
        array(i) = sum(citmap(:));
    end
    result = array;
end
function result = count_connections(obj)
% Method returning number of adjacent points in different regions.
% result = count_connections(obj)
% Return value is a row vector representing number of adjacent points for
% each region.
    counts = zeros(1,obj.regiN);
    for i = 1:obj.data.count
        neighbors = obj.data.neighborsOfPoint(i);
        neighbors = obj.regiList.var(neighbors,:);
        counts = counts + sum(neighbors & ~obj.regiList.var(i,:),1);
    end
    
    result = counts;
end
function result = roundness_of_regions(obj)
% Method regurning roundness of regions.
% result = roundness_of_regions(obj)
% Return value is a row vector with each element representing 
% 4 * pi * a / p^2 where p is a perimeter length of a region and a is
% an area of the region.
% Perimeter length is given by count_connections(), and so does not count
% boundary between regions and outside.
    areas = obj.size_regions;
    rims = obj.perimeter_length;
    result = 4 * pi * areas ./ (rims .^ 2);
end

function result = silhouette(obj)
% Method returning silhouette value of each point.
% result = silhouette(obj)
% Return value is a matrix with rows representing points and columns
% representing regions.
    silh = zeros(obj.data.count,obj.regiN);
    
    for i = 1:obj.data.count
        p = obj.regiList.var(i,:);
        if sum(obj.regiList.var(:,p)) < 2
            continue;
        end
        
        array = obj.meter.distanceArray(i);
        d = sum(array .* obj.regiList.var,1);
        a = d(p) / (sum(obj.regiList.var(:,p)) - 1);
        b = min(d(~p) ./ sum(obj.regiList.var(:,~p)));
        silh(i,p) = (b - a) / max([a,b]);
    end
    
    result = silh;
end

end
end
