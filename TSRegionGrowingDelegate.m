% TSRegionGrowingDelegate: objective class of delegate from
% run_region_growing().
% Written by Satoshi Yamashita.

classdef TSRegionGrowingDelegate < SYObject
properties
    hint = nan;
    loopCount = [];
end

methods
function dest = copy(obj,dest)
    if nargin < 2
        dest = TSRegionGrowingDelegate;
    end
    copy@SYObject(obj,dest);
    
    dest.hint = obj.hint;
    dest.loopCount = obj.loopCount;
end

function logLoopCount(obj,count)
    obj.loopCount = count;
end

end
end
