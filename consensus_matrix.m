
% Consensus matrix.
% Written by Satoshi Yamashita.

% This function returns a consensus matrix for given regions record.

function result = consensus_matrix(regiRec,hint)
% result = sampN x sampN matrix;

% Initialization.
regiRec = double(regiRec);

sampN = size(regiRec,1);
regiN = size(regiRec,2);
inseN = size(regiRec,3);

cm = zeros(sampN);

threshold = 0;
if ~isnan(hint.objectForKey('cm_threshold'))
    threshold = hint.objectForKey('cm_threshold');
end

for i = 1:inseN
    for j = 1:regiN
        cm = cm + regiRec(:,j,i) * regiRec(:,j,i)';
    end
end
cm = cm / inseN;

if threshold > 0
    cm = cm .* (cm > threshold);
end

result = cm;
end
