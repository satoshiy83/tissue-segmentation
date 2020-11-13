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
CM = SYData(consensus_matrix(regiRec,dict));
tnMatrix = SYData(zeros(3,length(20:2:80)));
tnMatrix.var(1,:) = (20:2:80) / 100;
tnMatrix.var(2,:) = nan;
tnMatrix.var(3,:) = 0;

%% Parameter search.
% Initial three points.
for t = [0.2,0.5,0.8]
    cm = CM.var .* (CM.var > t);
    partition = run_label_propagation(cm,hint);
    n = size(partition,2);
    index = tnMatrix.var(1,:) == t;
    tnMatrix.var(2,index) = n;
    tnMatrix.var(3,index) = 1;
end

stack = SYData(zeros(size(partition,1),regiN,iterN,'logical'));
counter = SYData(1);
loop = 0;
while loop < upthN
    i = find(tnMatrix.var(2,:) == regiN,1);
    if ~isempty(i)
        b = lf_try_label_propagation(tnMatrix,i,CM,hint, ...
            regiN,stack,counter,iterN);
        if b
            break;
        end
        continue;
    end
    
    i = find(tnMatrix.var(2,:) < regiN,1,'last');
    j = find(tnMatrix.var(2,:) > regiN,1,'first');
    if i + 1 == j || i == j + 1
        b = lf_try_label_propagation(tnMatrix,i,CM,hint, ...
            regiN,stack,counter,iterN);
        if b
            break;
        end
        b = lf_try_label_propagation(tnMatrix,j,CM,hint, ...
            regiN,stack,counter,iterN);
    elseif i > j
        i = round((i + j) / 2);
        b = lf_try_label_propagation(tnMatrix,i,CM,hint, ...
            regiN,stack,counter,iterN);
    else
        p = regiN - tnMatrix.var(2,i); q = tnMatrix.var(2,j) - regiN;
        i = round((i * q + j * p) / (p + q));
        b = lf_try_label_propagation(tnMatrix,i,CM,hint, ...
            regiN,stack,counter,iterN);
    end
    if b
        break;
    end
    
    loop = loop + 1;
end

result = stack.var;
end

function result = ...
    lf_try_label_propagation(tnMatrix,index,CM,hint, ...
        regiN,stack,counter,iterN)
cm = CM.var .* (CM.var > tnMatrix.var(1,index));
partition = run_label_propagation(cm,hint);
n = size(partition,2);
if isnan(tnMatrix.var(2,index))
    tnMatrix.var(2,index) = n;
else
    tnMatrix.var(2,index) = ...
        (tnMatrix.var(2,index) * tnMatrix.var(3,index) + n) / ...
        (tnMatrix.var(3,index) + 1);
end
tnMatrix.var(3,index) = tnMatrix.var(3,index) + 1;

if n == regiN
    stack.var(:,:,counter.var) = lf_sort_columns(partition);
    counter.var = counter.var + 1;
    if counter.var > iterN
        result = true;
        return;
    end
end
result = false;
end
function result = lf_sort_columns(partition)
array = zeros(1,size(partition,2));
for i = 1:size(partition,2)
    array(i) = find(partition(:,i),1);
end
[~,indices] = sort(array);

result = partition(:,indices);
end
