% Lable propagation clustering with specified number of clusters.
% Written by Satoshi Yamashita.

% This function implements label propagation algorithm and returns given
% number of clusters.

function result = run_cm_threshloding_lp(regiRec,hint)
% Function implementing label propagation algorithm.
% result = run_cm_threshloding_lp(regiRec,hint)
% Argument regiRec is a 3D matrix whose rows and columns represent data
% points and clusters, and slices represents different clustering results.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   lp_regiN: (double) number of clusters.
%   lp_upthN: (double) maximum number to update threshold.
%   lp_upcmN: (double) maximum number of label propagation feedback loop.
%   lp_iterN: (double) number of iteration in label propagation feedback
%   loop.
%   lp_uplaN: (double) maximum number of label update in 
%   run_label_propagation().
% Return value is a 3D logical matrix whose rows and columns represent data
% points and clusters and slices represent iterated trials.
% Each slice in the returned value is sorted.

%% Initialization.
regiN = hint.objectForKey('lp_regiN');
upthN = hint.objectForKey('lp_upthN');
upcmN = hint.objectForKey('lp_upcmN');
iterN = hint.objectForKey('lp_iterN');

if any(isnan([regiN,upthN,upcmN,iterN]))
    disp('Some hint is missing.');
    return;
end

dict = SYDictionary;
dict.setObjectForKey('cm_threshold',0.0);
CM = consensus_matrix(regiRec,dict);
array = (20:2:80) / 100;
array(2,:) = nan;
array(3,:) = 0;


%% Parameter search.
% Initial two points.
cm = CM .* (CM > 0.5);
partition = run_label_propagation(cm,hint);
n = size(partition,2);
if n == regiN
    result = partition;
    return;
end
array(2,array(1,:) == 0.5) = n;
array(3,array(1,:) == 0.5) = 1;
cm = CM .* (CM > 0.2);
partition = run_label_propagation(cm,hint);
n = size(partition,2);
if n >= regiN
    result = partition;
    return;
end
array(2,array(1,:) == 0.2) = n;
array(3,array(1,:) == 0.2) = 1;
cm = CM .* (CM > 0.8);
partition = run_label_propagation(cm,hint);
n = size(partition,2);
if n <= regiN
    result = partition;
    return;
end
array(2,array(1,:) == 0.8) = n;
array(3,array(1,:) == 0.8) = 1;

loop = 0;
while loop < upthN
    i = find(array(2,:) < regiN,1,'last');
    j = find(array(2,:) > regiN,1,'first');
    if i > j || i + 1 == j
        cm = CM .* (CM > array(1,i));
        partition = run_label_propagation(cm,hint);
        n = size(partition,2);
        if n == regiN
            break;
        end
        array(2,i) = (array(2,i) * array(3,i) + n) / (array(3,i) + 1);
        array(3,i) = array(3,i) + 1;
        cm = CM .* (CM > array(1,j));
        partition = run_label_propagation(cm,hint);
        n = size(partition,2);
        if n == regiN
            break;
        end
        array(2,j) = (array(2,j) * array(3,j) + n) / (array(3,j) + 1);
        array(3,j) = array(3,j) + 1;
    else
        p = regiN - array(2,i); q = array(2,j) - regiN;
        k = round((i * q + j * p) / (p + q));
        cm = CM .* (CM > array(1,k));
        partition = run_label_propagation(cm,hint);
        n = size(partition,2);
        if n == regiN
            break;
        end
        if isnan(array(2,k))
            array(2,k) = n;
        else
            array(2,k) = (array(2,k) * array(3,k) + n) / (array(3,k) + 1);
        end
        array(3,k) = array(3,k) + 1;
    end
    
    loop = loop + 1;
end

regiRec = zeros(size(partition,1),regiN,iterN,'logical');
for i = 1:iterN
    partition = run_label_propagation(cm,hint);
    partition = lf_sort_columns(partition);
    regiRec(:,1:size(partition,2),i) = partition;
end

result = regiRec;
end

function result = lf_sort_columns(partition)
array = zeros(1,size(partition,2));
for i = 1:size(partition,2)
    array(i) = find(partition(:,i),1);
end
[~,indices] = sort(array);

result = partition(:,indices);
end
