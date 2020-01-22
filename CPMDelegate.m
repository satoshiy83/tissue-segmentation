% Delegate class for run_cellular_Potts_model().
% Written by Satoshi Yamashita.

classdef CPMDelegate < SYObject
properties
    record = [];
    hint = nan;
end

methods
function dest = copy(obj,dest)
    if nargin < 2
        dest = CPMDelegate;
    end
    copy@SYObject(obj,dest);
    
    dest.record = obj.record;
    dest.hint = obj.hint;
end

end
end
