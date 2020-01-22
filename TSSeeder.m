% TSSeeder: objective class to make seeds.
% Written by Satoshi Yamashita.

classdef TSSeeder < TSRegionAnalyser
properties
%     data = nan;
    hint = nan;
    
    seedList = [];
%     regiList = [];
%     meter = nan;
    
%     regiN = [];
    collector = nan;
end

methods
function obj = TSSeeder(dataMap,hint)
% Tissue segmentation class to choose regions seeds.
% obj = TSSeeder(data,hint)
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   collector: (str) name of method collecting new seeds.
%   regiN: (double) number of regions.
    if nargin < 2
        return;
    end
    
    obj.initWithData(dataMap,hint);
end
function obj = initWithData(obj,dataMAp,hint)
% Initialization method with TSDataMap and SYDictionary instances.
% obj = initWithData(obj,data,hint)
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   collector: (str) name of method collecting new seeds.
%   regiN: (double) number of regions.
    initWithData@TSRegionAnalyser(obj,dataMAp);
    
    obj.hint = hint;
    
    str = 'centroids';
    if ~hint.isNanForKey('collector')
        str = hint.objectForKey('collector');
    end
    obj.collector = str2func(str);
    if ~hint.isNanForKey('rg_regiN') && isempty(obj.regiN)
        obj.regiN = hint.objectForKey('rg_regiN');
    end
end

function dest = copy(obj,dest)
% Method to make a copy.
% dest = copy(obj,dest)
    if nargin < 2
        dest = TSSeeder;
    end
    copy@RegionAnalyser(obj,dest);
    
    dest.seedList = obj.seedList;
    dest.collector = obj.collector;
    dest.regiN = obj.regiN;
end

function result = initialSeeds(obj)
% Method returning randomly chosen points.
% result = initialSeeds(obj)
% Return value is a row vector of points indices.
    array = randperm(obj.data.count);
    result = array(1:obj.regiN);
end
function result = initialMeans(obj)
% Method returning values of seeds.
% result = initialMeans(obj)
% Return value is a matrix with row vectors representing means.
    result = obj.data.dataList(obj.seedList.var,3:end);
end

function result = newSeeds(obj)
% Method returning new seeds.
% result = newSeeds(obj)
% Return value is a row vector of points indices.
    result = obj.collector(obj);
end
function result = newMeans(obj)
% Method returning new regions means.
% result = newMeans(obj)
% Return value is a matrix with row vectors representing means.
    result = obj.average_regions;
end

function result = centroids(obj)
% Method returning new seeds from regions centroids.
% result = centroids(obj)
    list = zeros(obj.regiN);
    for i = 1:obj.regiN
        points = obj.data.dataList(obj.regiList.var(:,i),1:2);
        aver = mean(points,1);
        d = sum((points - aver) .^ 2,2);
        index = find(d == min(d));
        indices = find(obj.regiList.var(:,i));
        if length(index) > 1
            points = obj.data.dataList(obj.regiList.var(:,i),3:end);
            aver = mean(points,1);
            d = obj.meter.distanceFromPoint(aver,indices(index));
            jndex = find(d == min(d));
            index = index(jndex(1));
        end
        list(i) = indices(index);
    end
    result = list;
end
function result = centroidsFirst(obj)
% Method returning new seeds from regions centroids.
% result = centroidsFirst(obj)
    list = zeros(obj.regiN);
    for i = 1:obj.regiN
        points = obj.data.dataList(obj.regiList.var(:,i),1:2);
        aver = mean(points,1);
        d = sum((points - aver) .^ 2,2);
        index = find(d == min(d));
        indices = find(obj.regiList.var(:,i));
        list(i) = indices(index(1));
    end
    result = list;
end

function result = initialQueue(obj)
% Method returning a queue to expand regions.
% result = initialQueue(obj)
% Return value is a (double) [index of point to expand cluster, index of
% point from which cluster expands, cluster id number, distance between
% point and mean, flag not used anymore].
    neigN = 8;
    if ~isnan(obj.hint.objectForKey('neigN'))
        neigN = obj.hint.objectForKey('neigN');
    end
    
    queue = zeros(obj.regiN * neigN,5);
    indices = 1:neigN;
    for i = 1:obj.regiN
        queue(indices,1) = obj.data.neigList(obj.seedList.var(i),:);
        queue(indices,2) = obj.seedList.var(i);
        queue(indices,3) = i;
        queue(indices,5) = true;
        
        indices = indices + neigN;
    end
    queue = queue(queue(:,1) > 0,:);
    for i = 1:size(queue,1)
        queue(i,4) = obj.meter.measure(queue(i,1:3));
    end
    queue = sortrows(queue,4);
    
    result = queue;
end

end
end
