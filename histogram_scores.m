% Function to convert an array of scalar values to a histogram.
% Written by Satoshi Yamashita.

function result = histogram_scores(scorArray,hint)
% Function to convert an array of scalar values to a histogram.
% result = histogram_scores(scorArray,hint)
% Argument scorArray is an array of scalar value.
% Argument hint is an SYDictionary insctance containing parameters.
% Return value is a matrix whose second row represents steps of value,
% third row represents number of elements in the array that is equal to or
% larger than a value in an upper cell and smaller than a right upper cell,
% and first row shows range and values dividing 5%, 0.5%, 0.05%, ..., and
% the rest of the values in scorArray.
% 
% Parameters:
%   hg_stepN: number of steps.
%   hg_sortD: sorting order to divide scorArray, either 'ascend' or
%   'descend'.
%   hg_range: range of histogram. If nan, the minimum and maximum values
%   are used to make the range.

% Initialization.
steps = 100;
if ~isnan(hint.objectForKey('hg_stepN'))
    steps = hint.objectForKey('hg_stepN');
end
sortD = 'ascend';
if ~isnan(hint.objectForKey('hg_sortD'))
    sortD = hint.objectForKey('hg_sortD');
end
if ~isnan(hint.objectForKey('hg_range'))
    range = hint.objectForKey('hg_range');
else
    range = [min(scorArray),max(scorArray)];
end
interval = (range(2) - range(1)) / (steps - 1);
hist = zeros(3,steps);

% Converting array.
hist(1,1:2) = range;
hist(2,1) = range(1);
for i = 2:steps
    hist(2,i) = hist(2,i - 1) + interval;
end
for i = 1:length(scorArray)
    index = floor((scorArray(i) - range(1)) / interval) + 1;
    index = min([index,steps]);
    hist(3,index) = hist(3,index) + 1;
end

scores = sort(scorArray,sortD);
divider = 20;
index = floor(length(scores) / divider);
p_values = [];
while index > 0
    p_values = cat(2,p_values,scores(index));
    divider = divider * 10;
    index = floor(length(scores) / divider);
end
hist(1,3:2 + length(p_values)) = p_values;

result = hist;
end
