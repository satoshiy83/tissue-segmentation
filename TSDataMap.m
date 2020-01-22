% TSDataMap: objective class containing tissue data to be segmented.
% Written by Satoshi Yamashita.

classdef TSDataMap < SYObject
properties
    data = nan;
    property_depth = 0;
    time_length = 0;
    hint = nan;
    
    dataList = [];
    neigList = [];
end

methods
function obj = TSDataMap(data,hint)
% Tissue segmentation data map class.
% obj = TSDataMap(data,hint)
% Argument data is a bitmap data.
% Argument hint is an SYDictionary instance containing parameters.
%
% Parameters:
%   property_depth: number of components per pixel at single slice.
%   time_length: number of time points.
%   neigN: number of neighboring pixels, 4 or 8.

    % Exceptional handlings.
    if nargin < 1
        return;
    elseif nargin == 1
        hint = SYDictionary;
    end
    
    obj.initWithData(data,hint);
end
function obj = initWithData(obj,data,hint)
% Initialization method with bitmap data.
% obj = initWithData(obj,data,hint)
% Argument data is a bitmap data.
% Argument hint is an SYDictionary instance containing parameters.
%
% Parameters:
%   property_depth: number of components per pixel at single slice.
%   time_length: number of time points.
%   neigN: number of neighboring pixels, 4 or 8.

    % Default initialization.
    if size(data,4) > 1
        data = ss_fold_timeSeries(data);
    end
    obj.data = data;
    obj.hint = hint;
    if ~isnan(hint.objectForKey('property_depth'))
        obj.property_depth = hint.objectForKey('property_depth');
    else
        obj.property_depth = size(data,3);
    end
    if ~isnan(hint.objectForKey('time_length'))
        obj.time_length = hint.objectForKey('time_length');
    else
        obj.time_length = floor(size(data,3) / obj.property_depth);
    end
    if size(data,3) ~= obj.property_depth * obj.time_length
        disp('Bit-depth mismatch!');
    end
    
    obj.dataList = lf_list_data(obj);
    obj.neigList = lf_list_neighbors(obj);
end
function obj = initWithDataList(obj,dataList,frame,hint)
% Initialization method with listed data.
% obj = initWithDataList(obj,dataList,frame,hint)
% Argument dataList is a list of data (double) [x,y,value...].
% Argument frame specifies bitmap size double[w,h].
% Argument hint is an SYDictionary instance containing parameters.
%
% Parameters:
%   property_depth: number of components per pixel at single slice.
%   time_length: number of time points.
%   neigN: number of neighboring pixels, 4 or 8.

    obj.dataList = dataList;
    obj.hint = hint;
    obj.data = lf_delist_dataList(obj,frame);
    
    if ~isnan(hint.objectForKey('property_depth'))
        obj.property_depth = hint.objectForKey('property_depth');
    else
        obj.property_depth = size(obj.data,3);
    end
    if ~isnan(hint.objectForKey('time_length'))
        obj.time_length = hint.objectForKey('time_length');
    else
        obj.time_length = floor(size(obj.data,3) / obj.property_depth);
    end
    if size(obj.data,3) ~= obj.property_depth * obj.time_length
        disp('Bit-depth mismatch!');
    end
    
    obj.neigList = lf_list_neighbors(obj);
end

function dest = copy(obj,dest)
% Method to make a copy.
% dest = copy(obj,dest)
    if nargin < 2
        dest = TSDataMap;
    end
    copy@SYObject(obj,dest);

    dest.data = obj.data;
    dest.property_depth = obj.property_depth;
    dest.time_length = obj.time_length;
    
    dest.dataList = obj.dataList;
    dest.neigList = obj.neigList;
end

function result = frameSize(obj)
% Method returning frame size.
% result = frameSize(obj)
    frame = size(obj.data);
    result = frame(1:2);
end

function pbj = sliceAtTime(obj,t,flag_timeMask)
% Method returning TSDataMap instance of given time point.
% pbj = sliceAtTime(obj,t,flag_timeMask)
% Argument t specifies a time point.
% Argument flag_timeMask is a boolean indicating masking pixels having nan
% value at any time point.
    if t < 1 || t > obj.time_length
        pbj = nan;
        return;
    end
    
    indices = (1:obj.property_depth) + obj.property_depth * (t - 1);
    slice = obj.data(:,:,indices);
    if flag_timeMask
        mask = ones(size(obj.data(:,:,1)));
        mask(any(isnan(obj.data),3)) = nan;
        slice = slice .* mask;
    end
    dict = obj.hint.copy;
    dict.setObjectForKey('time_length',1);
    
    pbj = TSDataMap(slice,dict);
end

function result = count(obj)
% Method counting samples in the data.
% result = count(obj)
    result = size(obj.dataList,1);
end
function result = indexOfPoint(obj,point)
% Method to return index of point (x,y) in the data list.
% result = indexOfPoint(obj,point)
    result = find(all(obj.dataList(:,1:2) == point,2));
end

function result = dataAtPoint(obj,point)
% Method returning data at given pixel.
% result = dataAtPoint(obj,point)
% Argument point is either an index, point, or point-time.
    result = nan;
    if length(point) == 1
        result = obj.dataList(point,3:end);
    elseif length(point) == 2
        index = all(obj.dataList(:,1:2) == point,2);
        result = obj.dataList(index,3:end);
%         result = permute(obj.data(point(2),point(1),:),[1,3,2]);
    elseif length(point) == 3
        indices = 1:obj.property_depth;
        indices = indices + obj.property_depth * (array(3) - 1);
        index = all(obj.dataList(:,1:2) == point(1:2),2);
        result = obj.dataList(index,indices);
%         result = permute(obj.data(point(2),point(1),indices),[1,3,2]);
    end
end

function result = unfoldData(obj)
% Method to unfold 3D matrix to 4D.
% result = unfoldData(obj)
    
    frameSize = obj.frameSize;
    matrix = zeros(frameSize(1),frameSize(2), ...
        obj.property_depth,obj.time_length);
    range = 1:obj.property_depth;
    for t = 1:obj.time_length
        matrix(:,:,:,t) = obj.data(:,:,range);
        range = range + obj.property_depth;
    end
    
    result = matrix;
end

function result = neighborsOfPoint(obj,point)
% Method returning points neighboring to a given point.
% result = neighborsOfPoint(obj,point)
    if length(point) == 1
        index = point;
    elseif length(point) == 2
        index = find(all(obj.dataList(:,1:2) == point,2));
    else
        result = nan;
        return;
    end
    
    array = obj.neigList(index,:);
    result = array(array > 0);
end

end
end

% Local functions.
function result = lf_list_data(obj)
data_size = size(obj.data);
if length(data_size) < 3
    data_size = [data_size,1,1];
end
ny = data_size(1); nx = data_size(2);
nBoxes = nx * ny;

list = zeros(nBoxes,2 + data_size(3));
i = 0;
for b = 1:nBoxes
    [ky,kx] = ind2sub([ny,nx],b);
    if any(isnan(obj.data(ky,kx,:)))
        continue;
    end
    i = i + 1;
    list(i,:) = [kx,ky,permute(obj.data(ky,kx,:),[1,3,2])];
end

result = list(1:i,:);
end
function result = lf_delist_dataList(obj,frame)
data = nan(frame(3),frame(4),size(obj.dataList,2) - 2);
for i = 1:size(obj.dataList,1)
    data(obj.dataList(i,2),obj.dataList(i,1),:) = obj.dataList(i,3:end);
end

result = data;
end
function result = lf_list_neighbors(obj)
neigN = 8;
if ~isnan(obj.hint.objectForKey('neigN'))
    neigN = obj.hint.objectForKey('neigN');
end
list = zeros(obj.count,neigN);

if neigN == 4
    for i = 1:obj.count
        x = obj.dataList(i,1); y = obj.dataList(i,2);
        x = x - 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,1) = index;
        end
        x = x + 1;
        y = y - 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,2) = index;
        end
        y = y + 1;
        x = x + 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,3) = index;
        end
        x = x - 1;
        y = y + 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,4) = index;
        end
    end
elseif neigN == 8
    for i = 1:obj.count
        x = obj.dataList(i,1); y = obj.dataList(i,2);
        x = x - 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,1) = index;
        end
            y = y - 1;
            index = find(all(obj.dataList(:,1:2) == [x,y],2));
            if ~isempty(index)
                list(i,2) = index;
            end
            y = y + 1;
            y = y + 1;
            index = find(all(obj.dataList(:,1:2) == [x,y],2));
            if ~isempty(index)
                list(i,3) = index;
            end
            y = y - 1;
        x = x + 1;
        y = y - 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,4) = index;
        end
        y = y + 1;
        x = x + 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,5) = index;
        end
            y = y - 1;
            index = find(all(obj.dataList(:,1:2) == [x,y],2));
            if ~isempty(index)
                list(i,6) = index;
            end
            y = y + 1;
            y = y + 1;
            index = find(all(obj.dataList(:,1:2) == [x,y],2));
            if ~isempty(index)
                list(i,7) = index;
            end
            y = y - 1;
        x = x - 1;
        y = y + 1;
        index = find(all(obj.dataList(:,1:2) == [x,y],2));
        if ~isempty(index)
            list(i,8) = index;
        end
    end
end

result = list;
end
