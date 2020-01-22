% Label propagation clustering.
% Written by Satoshi Yamashita.

% This function implements community detection algorithm using label
% propagation.

function result = run_label_propagation(cm,hint)
% Function implementing label propagation algorithm.
% result = run_label_propagation(cm,hint)
% Argument cm is an adjacency matrix of weighted graph in which light edges
% were already removed.
% Argument hint is an SYDictionary instance holding parameterse.
% Parameters:
%   lp_uplaN: (double) maximum number of label update.
% Return value is a logical matrix whose rows and columns represent data
% points and clusters.

%% Initialization.
sampN = size(cm,1);
I = eye(sampN) == 0;
cm = cm .* I;
labels = 1:sampN;
labels = labels';

uplaN = hint.objectForKey('lp_uplaN');

regiRec = [];

%% Main label propagation loop.
while ~(uplaN < 1)
    flag = true;
    for i = randperm(sampN)
        [~,indices] = max(sum(cm(:,i) .* (labels == 1:sampN),1));
        if labels(i) ~= indices(1)
            flag = false;
        end
        labels(i) = indices(1);
    end
    
    if flag
        break;
    end
    
    uplaN = uplaN - 1;
end

for i = 1:sampN
    m = max(labels);
    region = labels == m;
    regiRec = cat(2,regiRec,region);
    
    labels = labels .* (region == 0);
    if sum(labels) == 0
        break;
    end
end

result = logical(regiRec);
end
