% TSMeter: objective class to measure distance in property space.
% Written by Satoshi Yamashita.

% TSMeter calculate distance between points in TSDataMap.
% For subclasses of TSMeter, metric methods should follow the rules below.
% For a metric named "base", TSMeter subclass should implement methods
% function result = base(obj,point,qoint),
% function result = base_adjacent(obj,array),
% function result = base_array(obj,point),
% and can implement method with any extension (as "ext" in below)
% function result = base_ext(obj,array).
% The arguments point and qoint are point values.
% The argument array is a double array [index of point to expand cluster, 
% index of point from which cluster expands, cluster id number].


classdef TSMeter < SYObject
properties
    data = nan; % (TSDataMap) given at initialization.
    hint = nan; % (SYDictionary) given at initialization.
    
    seedList = []; % (SYData, shared) given after init.
    meanList = []; % (SYData, shared) given after init.
    regiList = []; % (SYData, shared) given after init.
    
    metric = nan; % (function_handle) given in the hint.
    metric_base = nan; % (function_handle) given in the hint.
    metric_adjacent = nan; % (function_handle) given in the hint.
    metric_array = nan; % (function_handle) given in the hint.
end

methods
function obj = TSMeter(data,hint)
% Tissue segmentation class defining metric.
% obj = TSMeter(data,hint)
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% 
% Parameters:
%   metric_base: part of method name, defining distance between two points.
%   metric_exte: part of method name, defining distance between point and
%   cluster.

    if nargin < 2
        return;
    end
    
    obj.initWithData(data,hint);
end
function obj = initWithData(obj,data,hint)
% Initialization method with TSDataMap and hint.
% obj = initWithData(obj,data,hint)
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% 
% Parameters:
%   metric_base: part of method name, defining distance between two points.
%   metric_exte: part of method name, defining distance between point and
%   cluster.

    obj.data = data;
    obj.hint = hint;
    
    str_base = 'Euclidean';
    if ~isnan(hint.objectForKey('metric_base'))
        str_base = hint.objectForKey('metric_base');
    end
    str_exte = 'mean';
    if ~isnan(hint.objectForKey('metric_exte'))
        str_exte = hint.objectForKey('metric_exte');
    end
    obj.metric = str2func([str_base,'_',str_exte]);
    obj.metric_base = str2func(str_base);
    obj.metric_adjacent = str2func([str_base,'_adjacent']);
    obj.metric_array = str2func([str_base,'_array']);
end

function dest = copy(obj,dest)
% Method to make a copy.
% dest = copy(obj,dest)
    if nargin < 2
        dest = TSMeter;
    end
    copy@SYObject(obj,dest);
    
    dest.data = obj.data;
    dest.hint = obj.hint;
    dest.seedList = obj.seedList;
    dest.meanList = obj.meanList;
    dest.regiList = obj.regiList;
    dest.metric = obj.metric;
    dest.metric_adjacent = obj.metric_adjacent;
    dest.metric_array = obj.metric_array;
end

function result = measure(obj,array)
% Method returning distance between a point and cluster.
% result = measure(obj,array)
% Argument array is (double) [index of point to expand cluster, index of
% point from which cluster expands, cluster id number].
    result = obj.metric(obj,array);
end
function result = distanceBetweenPoints(obj,point,qoint)
% Method returning distance between two points.
% result = distanceBetweenPoints(obj,point,qoint)
% Arguments point and qoint are indices of points.
    result = obj.metric_base(obj,point,qoint);
end
function result = distanceMatrix(obj,array,brray)
% Method returning distances between points in two array.
% result = distanceMatrix(obj,array,brray)
% Arguments (double[n]) array and brray contain indices of points.
    table = zeros(length(array),length(brray));
    for i = 1:length(array)
        for j = 1:length(brray)
            point = obj.data.dataList(array(i),3:end);
            qoint = obj.data.dataList(brray(j),3:end);
            table(i,j) = obj.metric_base(point,qoint);
        end
    end
    result = table;
end
function result = distanceArray(obj,point)
% Method returning distance between a given point and other points.
% result = distanceArray(obj,point)
% Argument point is an index of the point.
    if isnan(obj.hint.objectForKey('TSMeterPointConversionSkip')) || ...
            ~obj.hint.objectForKey('TSMeterPointConversionSkip')
        if length(point) < 4 && all(point == floor(point))
            point = obj.data.dataAtPoint(point);
        end
    end
    
    result = obj.metric_array(obj,point);
end
function result = distanceFromPoint(obj,point,array)
% Method returning distance between a given point and points in an array.
% result = distanceFromPoint(obj,point,array)
% Argument point is a point value.
% Argument (double[n]) is an array of indices.
    list = zeros(length(array),1);
    for i = 1:length(array)
        qoint = obj.data.dataList(array(i),3:end);
        list(i) = obj.metric_base(obj,point,qoint);
    end
    result = list;
end

function result = absolute(~,point,qoint)
% Metric of L-1 norm.
% result = absolute(~,point,qoint)
    result = abs(point - qoint);
end
function result = absolute_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    result = abs(a - b);
end
function result = absolute_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    result = abs(a - b);
end
function result = absolute_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    result = abs(a - b);
end
function result = absolute_array(obj,point)
    a = point;
    d = obj.data.dataList(:,3:end) - a;
    result = abs(d);
end

function result = time_integral_absolute(~,point,qoint)
% Metric of L-1 norm.
% time_integral_absolute(obj,point,qoint)
    result = sum(abs(point - qoint));
end
function result = time_integral_absolute_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    result = sum(abs(a - b));
end
function result = time_integral_absolute_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    result = sum(abs(a - b));
end
function result = time_integral_absolute_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    result = sum(abs(a - b));
end
function result = time_integral_absolute_array(obj,point)
    a = point;
    d = obj.data.dataList(:,3:end) - a;
    result = sum(abs(d),2);
end

function result = Euclidean(~,point,qoint)
% Metric of L-2 norm.
% result = Euclidean(obj,point,qoint)
    result = sqrt(sum((point - qoint) .^ 2));
end
function result = Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    result = sqrt(sum((a - b) .^ 2));
end
function result = Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    result = sqrt(sum((a - b) .^ 2));
end
function result = Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    result = sqrt(sum((a - b) .^ 2));
end
function result = Euclidean_array(obj,point)
    a = point;
    d = obj.data.dataList(:,3:end) - a;
    result = sqrt(sum(d .^ 2,2));
end

function result = time_integral_Euclidean(obj,point,qoint)
% Metric of L-2 norm summed through time.
% result = time_integral_Euclidean(obj,point,qoint)
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = (point - qoint) .^ 2;
    result = sum(sqrt(sum(table,1)));
end
function result = time_integral_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = (a - b) .^ 2;
    result = sum(sqrt(sum(table,1)));
end
function result = time_integral_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = (a - b) .^ 2;
    result = sum(sqrt(sum(table,1)));
end
function result = time_integral_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = (a - b) .^ 2;
    result = sum(sqrt(sum(table,1)));
end
function result = time_integral_Euclidean_array(obj,point)
    a = point;
    d = obj.data.dataList(:,3:end) - a;
    
    d = d .^ 2;
    s = zeros(obj.data.count,1);
    indices = 1:obj.data.property_depth;
    for t = 1:obj.data.time_length
        s = s + sqrt(sum(d(:,indices),2));
        indices = indices + obj.data.property_depth;
    end
    result = s;
end

function result = tensor_Euclidean(~,tensor,uensor)
% Metric of tensor Euclidean.
% result = tensor_Euclidean(obj,tensor,uensor)
    point = [(tensor(1) - tensor(4)) / 2,tensor(2)];
    qoint = [(uensor(1) - uensor(4)) / 2,uensor(2)];
    result = sqrt(sum((point - qoint) .^ 2));
end
function result = tensor_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    a = [(a(1) - a(4)) / 2,a(2)];
    b = [(b(1) - b(4)) / 2,b(2)];
    result = sqrt(sum((a - b) .^ 2));
end
function result = tensor_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    a = [(a(1) - a(4)) / 2,a(2)];
    b = [(b(1) - b(4)) / 2,b(2)];
    result = sqrt(sum((a - b) .^ 2));
end
function result = tensor_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    a = [(a(1) - a(4)) / 2,a(2)];
    b = [(b(1) - b(4)) / 2,b(2)];
    result = sqrt(sum((a - b) .^ 2));
end
function result = tensor_Euclidean_array(obj,point)
    a = point;
    a = [(a(1) - a(4)) / 2,a(2)];
    d = obj.data.dataList(:,3:end);
    d = [(d(:,1) - d(:,4)) / 2,d(:,2)] - a;
    
    result = sqrt(sum(d .^ 2,2));
end

function result = time_integral_tensor_Euclidean(obj,point,qoint)
% Metric of tensor Euclidean summed through time.
% result = time_integral_tensor_Euclidean(obj,point,qoint)
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = point;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = qoint;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_tensor_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = b;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_tensor_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = b;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_tensor_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = b;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_tensor_Euclidean_array(obj,point)
    a = (point(1:4:end) - point(4:4:end)) / 2;
    b = point(2:4:end);
    c = (obj.data.dataList(:,3:4:end) - obj.data.dataList(:,6:4:end)) / 2;
    d = obj.data.dataList(:,4:4:end);
    
    result = sum(sqrt((a - c) .^ 2 + (b - d) .^ 2),2);
end

function result = stretched_linear_Euclidean(~,tensor,uensor)
% Metric of L-2 norm.
% result = stretched_linear_Euclidean(~,tensor,uensor)
    result = sqrt(sum((tensor - uensor) .^ 2));
end
function result = stretched_linear_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    c = a(1:4);
    d = b(1:4) + b(5) * a(5:8);
    result = sqrt(sum((c - d) .^ 2));
end
function result = stretched_linear_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    result = sqrt(sum((a - b) .^ 2));
end
function result = stretched_linear_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    result = sqrt(sum((a - b) .^ 2));
end
function result = stretched_linear_Euclidean_array(obj,point)
    if length(point) == 5
        d = obj.data.dataList(:,3:end);
        d = point(1:4) + point(5) * d(:,5:8) - d(:,1:4);
    else
        a = point;
        d = obj.data.dataList(:,3:end) - a;
    end
    
    result = sqrt(sum(d .^ 2,2));
end

function result = stretched_linear_tensor_Euclidean(~,tensor,uensor)
% Metric of tensor Euclidean.
% result = stretched_linear_tensor_Euclidean(obj,tensor,uensor)
    point = [(tensor(1) - tensor(4)) / 2,tensor(2),(tensor(5) - tensor(8)) / 2,tensor(6)];
    qoint = [(uensor(1) - uensor(4)) / 2,uensor(2),(uensor(5) - uensor(8)) / 2,uensor(6)];
    result = sqrt(sum((point - qoint) .^ 2));
end
function result = stretched_linear_tensor_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    c = [(a(1) - a(4)) / 2,a(2)];
    d = [(b(1) - b(4)) / 2,b(2)] + b(5) * [(a(5) - a(8)) / 2,a(6)];
    result = sqrt(sum((c - d) .^ 2));
end
function result = stretched_linear_tensor_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    a = [(a(1) - a(4)) / 2,a(2),(a(5) - a(8)) / 2,a(6)];
    b = [(b(1) - b(4)) / 2,b(2),(b(5) - b(8)) / 2,b(6)];
    result = sqrt(sum((a - b) .^ 2));
end
function result = stretched_linear_tensor_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    a = [(a(1) - a(4)) / 2,a(2),(a(5) - a(8)) / 2,a(6)];
    b = [(b(1) - b(4)) / 2,b(2),(b(5) - b(8)) / 2,b(6)];
    result = sqrt(sum((a - b) .^ 2));
end
function result = stretched_linear_tensor_Euclidean_array(obj,point)
    if length(point) == 5
        d = obj.data.dataList(:,3:end);
        d = point(1:4) + point(5) * d(:,5:8) - d(:,1:4);
        d = [(d(:,1) - d(:,4)) / 2,d(:,2)];
    else
        a = point;
        a = [(a(1) - a(4)) / 2,a(2),(a(5) - a(8)) / 2,a(6)];
        d = obj.data.dataList(:,3:end);
        d = [(d(:,1) - d(:,4)) / 2,d(:,2),(d(:,5) - d(:,8)) / 2,d(:,6)] - a;
    end
    
    result = sqrt(sum(d .^ 2,2));
end

function result = time_integral_stretched_linear_tensor_Euclidean(obj,tensor,uensor)
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = tensor;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = uensor;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:); (table(5,:) - table(8,:)) / 2; table(6,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:); (uable(5,:) - uable(8,:)) / 2; uable(6,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_stretched_linear_tensor_Euclidean_mean(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.meanList.var(array(3),:);
    
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:)];
    uable = zeros(5,obj.data.time_length);
    uable(:) = b;
    uable(1:4,:) = uable(1:4,:) + uable(5,:) .* table(5:8,:);
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_stretched_linear_tensor_Euclidean_seed(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(obj.seedList.var(array(3)));
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = b;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:); (table(5,:) - table(8,:)) / 2; table(6,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:); (uable(5,:) - uable(8,:)) / 2; uable(6,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_stretched_linear_tensor_Euclidean_adjacent(obj,array)
    if array(1) == 0
        result = 0;
        return;
    end
    
    a = obj.data.dataAtPoint(array(1));
    b = obj.data.dataAtPoint(array(2));
    table = zeros(obj.data.property_depth,obj.data.time_length);
    table(:) = a;
    uable = zeros(obj.data.property_depth,obj.data.time_length);
    uable(:) = b;
    a = [(table(1,:) - table(4,:)) / 2; table(2,:); (table(5,:) - table(8,:)) / 2; table(6,:)];
    b = [(uable(1,:) - uable(4,:)) / 2; uable(2,:); (uable(5,:) - uable(8,:)) / 2; uable(6,:)];
    
    result = sum(sqrt(sum((a - b) .^ 2,1)));
end
function result = time_integral_stretched_linear_tensor_Euclidean_array(obj,point)
    a = (point(1:4:end) - point(4:4:end)) / 2;
    b = point(2:4:end);
    c = (obj.data.dataList(:,3:4:end) - obj.data.dataList(:,6:4:end)) / 2;
    d = obj.data.dataList(:,4:4:end);
    
    result = sum(sqrt((a - c) .^ 2 + (b - d) .^ 2),2);
end

end
end
