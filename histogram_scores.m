% Function to convert an array of scalar values to a histogram.
% Written by Satoshi Yamashita.

function result = histogram_scores(scorArray,hint)

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
