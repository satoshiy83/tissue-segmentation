% Function to give random segmentation for given times.
% Written by Satoshi Yamashita.

function result = run_nonsense_segmentation(dataMap,hint)
% Function to segment tissue randomly for given times.
% result = run_nonsense_segmentation(dataMap,hint)
% Argument dataMap is a TSDataMap instance.
% Argument hint iss an SYDictionary instance holding parameters.
% Return value is a matrix whose column vector shows the score of the
% segmentation.
% 
% Parameters:
%   ns_regiN: (int) number of segments.
%   ns_bootN: (int) number of trials.
%   ns_insoC: (int) minimum neighborhood size around initial seeds.
%   ns_fh_score: (char) name of function to score the segmentations.
% 
% Local functions:
%   result = lf_silhouette(analyser)
%       It returns total silhouette value.

% Initialization.
regiN = hint.objectForKey('ns_regiN');

bootN = 20000;
if ~isnan(hint.objectForKey('ns_bootN'))
    bootN = hint.objectForKey('ns_bootN');
end
insoC = 1;
if ~isnan(hint.objectForKey('ns_insoC'))
    insoC = hint.objectForKey('ns_insoC');
end
neigN = size(dataMap.neigList,2);
insoC = insoC * regiN * neigN;

fh_score = @lf_silhouette;
if ~isnan(hint.objectForKey('ns_fh_score'))
    fh_score = str2func(hint.objectForKey('ns_fh_score'));
end
meter = TSMeter(dataMap,hint);
analyser = TSRegionAnalyser(dataMap);
analyser.meter = meter;
partData = SYData;
meter.regiList = partData;
analyser.regiList = partData;

seeder = TSSeeder(dataMap,hint);
alloList = zeros(dataMap.count,1);
colu = fh_score(analyser);
scorTable = zeros(length(colu),bootN);

% Main loop of clustering.
t0 = datetime;
disp(['Main loop run from ',char(t0)]);
messanger = dialog('Position',[300 300 200 100],'WindowStyle','normal', ...
    'Name','k-mean processing');
msg1 = uicontrol(messanger,'Style','text','Position',[5 5 190 40], ...
    'FontSize',14);
for t = 1:bootN
    % Displaying estimated duration of calculation.
    if t == 1
        t1 = datetime;
        edc = (t1 - t0) * (bootN - 1);
        string = [char(t1),'; EDC: ',char(edc)];
        disp(string);
    elseif t == 10
        t1 = datetime;
        edc = (t1 - t0) * (bootN - 10) / 10;
        string = [char(t1),'; EDC: ',char(edc)];
        disp(string);
    end
    if ishandle(messanger)
        msg1.String = ['Iterating ',num2str(t),' / ',num2str(bootN)];
        drawnow;
    end

    % Initialization.
    seedList = seeder.initialSeeds;
    alloList(:) = 0;
    partition = zeros(dataMap.count,regiN,'logical');
    queue = zeros(regiN * neigN,4);
    indices = 1:neigN;
    for i = 1:regiN
        alloList(seedList(i)) = 1;
        partition(seedList(i),i) = true;
        
        queue(indices,1) = dataMap.neigList(seedList(i),:);
        queue(indices,2) = seedList(i);
        queue(indices,3) = i;
        queue(indices,4) = true;
        
        indices = indices + neigN;
    end
    count = insoC;
    
    % Random growing of regions.
    while size(queue,1) > 0
        if queue(1,1) > 0 && alloList(queue(1,1)) < 1
            partition(queue(1,1),queue(1,3)) = 1;
            alloList(queue(1,1)) = 1;
            array = dataMap.neighborsOfPoint(queue(1,1));
            push = zeros(length(array),4);
            j = 0;
            for i = 1:length(array)
                if alloList(array(i)) < 1
                    j = j + 1;
                    push(j,1) = array(1,i); % point to expand.
                    push(j,2) = queue(1,1); % point of expansion.
                    push(j,3) = queue(1,3); % region id.
                    push(j,4) = true;
                end
            end
            queue(1,:) = [];
            if j > 0
                queue = cat(1,queue,push(1:j,:));
                if count < 1
                    queue = queue(randperm(size(queue,1)),:);
                else
                    count = count - 1;
                end
            end
        else
            queue(1,:) = [];
        end
    end
    
    partData.var = partition;
    scorTable(:,t) = fh_score(analyser);
%     if data.time_length > 1
%         scorTable(:,t) = silhouette_atAllTimes(partition,data);
%     else
%         analyser.regiList = SYData(partition);
%         scorTable(1,t) = sum(sum(analyser.silhouette));
%     end
end

if ishandle(messanger)
    delete(messanger);
end

result = scorTable;
end

function result = lf_silhouette(analyser)
if isempty(analyser.regiList.var)
    result = 0;
    return;
end

result = sum(sum(analyser.silhouette));
end
function result = lf_silhouette_atAllTimes(analyser)
if isempty(analyser.regiList.var)
    result = zeros(2 + analyser.time_length,1);
    return;
end

result = silhouette_atAllTimes(analyser.regiList.var,analyser.data);
end
