% Iterated label propagation clustering.
% Written by Satoshi Yamashita.

% This function run label propagation first on a consensus matrix of given
% partitions, then resultant partitions are converted to a new consensus
% matrix and feeded to label propagation again. This process is iterated
% until it reaches a convergence.

function result = run_iterated_lp(regiRec,hint)
% Function implementing label propagation on a consensus matrix until it
% reaches convergence.
% result = run_iterated_lp(regiRec,hint)
% Argument regiRec is a 3D matrix whose rows and columns represent data
% points and clusters, and slices represents different clustering results.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   upseN: (double) maximum number to iterate label propagation.
%   upthN: (double) maximum number to update threshold for
%   run_cm_thresholding_lp().
%   iterN: (double) maximum number of label update for 
%   run_label_propagation().
% Return value is a logical matrix whose rows and columns represent data
% points and clusters.

%% Initialization.
regiN = size(regiRec,2);
iterN = size(regiRec,3);
upseN = hint.objectForKey('upseN');

hint = hint.copy;
hint.setObjectForKey('regiN',regiN);

%% Iterate label propagation.
segiRec = zeros(size(regiRec),'logical');
for t = 1:upseN
    segiRec(:) = false;
    for i = 1:iterN
        partition = run_cm_threshloding_lp(regiRec,hint);
        partition = lf_sort_columns(partition);
        segiRec(:,1:size(partition,2),i) = partition;
    end
    
    regiRec = segiRec == partition;
    if all(regiRec(:))
        break;
    end
    
    regiRec = segiRec;
end

disp('Iterated label propagations ',num2str(t),' times.');

result = partition;
end

function result = lf_sort_columns(partition)
array = zeros(1,size(partition,2));
for i = 1:size(partition,2)
    array(i) = find(partition(:,i),1);
end
[~,indices] = sort(array);

result = partition(:,indices);
end
