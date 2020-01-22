% Region growing segmentation.
% Written by Satoshi Yamashita.
% 
% This function run a segmentation by region growing algorithm in 2D bitmap
% of any kind of data-format, bit-depth.

function result = run_region_growing(dataMap,hint,delegate)
% Function implementing region growing algorithm.
% result = run_region_growing(data,hint,delegate)
% Argument data is a TSDataMap instance.
% Argument hint is a SYDictionary instance containing parameters.
% Argument delegate is a TSRegionGrowingDelegate instance.
% Return value is a double[n,k,t], where n is a size of sample, k is a
% number of partitions, and t is a number of trials.
%
% Parameters:
%   rg_regiN: number of regions.
%   rg_inseN: number of trials.
%   rg_upseN: maximum number of seeds update.
%   rg_cyclN: number to check periodic seeding.

%% Initialization.
if ~isa(dataMap,'TSDataMap')
    dataMap = TSDataMap(dataMap,hint);
end
regiN = hint.objectForKey('rg_regiN');
inseN = hint.objectForKey('rg_inseN');
upseN = hint.objectForKey('rg_upseN');
cyclN = hint.objectForKey('rg_cyclN');
meter = TSMeter(dataMap,hint);
seeder = TSSeeder(dataMap,hint);
if any(isnan([regiN,inseN,upseN,cyclN])) ...
        || isnan(dataMap) || isnan(meter) || isnan(seeder)
    disp('Some hints were missing!');
    result = nan;
    return;
end

alloList = zeros(dataMap.count,1);
seedList = SYData;
meanList = SYData;
regiList = SYData(zeros(dataMap.count,regiN,'logical'));
meter.seedList = seedList;
meter.meanList = meanList;
meter.regiList = regiList;
seeder.seedList = seedList;
seeder.regiList = regiList;
seeder.meter = meter;
delegate.hint = hint;

regiRec = [];

%% Main loop of clustering.
t0 = datetime;
disp(['Main loop run from ',char(t0)]);
messanger = dialog('Position',[300 300 200 100],'WindowStyle','normal', ...
    'Name','k-mean processing');
msg1 = uicontrol(messanger,'Style','text','Position',[5 5 190 40], ...
    'FontSize',14);
msg2 = uicontrol(messanger,'Style','text','Position',[5 55 190 40], ...
    'FontSize',14);
fig = figure;
frame = dataMap.frameSize;
for t = 1:inseN
    % Displaying estimated duration of calculation.
    if t == 1
        t1 = datetime;
        edc = (t1 - t0) * (inseN - 1);
        string = [char(t1),'; EDC: ',char(edc)];
        disp(string);
    elseif t == 10
        t1 = datetime;
        edc = (t1 - t0) * (inseN - 10) / 10;
        string = [char(t1),'; EDC: ',char(edc)];
        disp(string);
    end
    if ishandle(messanger)
        msg1.String = ['Iterating ',num2str(t),' / ',num2str(inseN)];
        drawnow;
    end
    
    % Initialization of seeds list and means list.
    seedList.var = seeder.initialSeeds;
    meanList.var = seeder.initialMeans;
    
    cyclC = cyclN;
    cyclMeanList = meanList.var;
    
    % Loop of updating seeds.
    for loop = 1:upseN
        if ishandle(messanger)
            msg2.String = ['Loop in ',num2str(loop),' / ',num2str(upseN)];
            drawnow;
        end
        
        % Initialization.
        regiList.var(:) = 0;
        alloList(:) = 0;
        for i = 1:regiN
            regiList.var(seedList.var(i),i) = 1;
            alloList(seedList.var(i)) = 1;
        end
        queue = seeder.initialQueue;
        
        % Expansion of regions.
        while size(queue,1) > 0
            if queue(1,1) > 0 && alloList(queue(1,1)) < 1
                regiList.var(queue(1,1),queue(1,3)) = 1;
                alloList(queue(1,1)) = 1;
                array = dataMap.neighborsOfPoint(queue(1,1));
                push = zeros(length(array),5);
                j = 0;
                for i = 1:length(array)
                    if alloList(array(i)) < 1
                        j = j + 1;
                        push(j,1) = array(1,i); % point to expand.
                        push(j,2) = queue(1,1); % point of expansion.
                        push(j,3) = queue(1,3); % region id.
                        push(j,4) = meter.measure(push(j,1:3));
                        push(j,5) = true;
                    end
                end
                queue(1,:) = [];
                if j > 0
                    queue = cat(1,queue,push(1:j,:));
                    queue = sortrows(queue,4);
                end
            else
                queue(1,:) = [];
            end
        end
        
        % Recording means to check convergence.
        lastMeanList = meanList.var;
        if cyclC > 0
            cyclC = cyclC - 1;
        elseif cyclC == 0
            cyclMeanList = cat(3,cyclMeanList,meanList.var);
            cyclC = cyclN;
        end
        % Updating seeds.
        seedList.var = seeder.newSeeds;
        meanList.var = seeder.newMeans;
        % Checking convergence.
        if isequal(lastMeanList,meanList.var) || ...
                any(all(all(cyclMeanList == meanList.var,1),2),3)
            break;
        end
        
        % drawing partition.
        if ishandle(fig)
            image = regiList2image(regiList.var,dataMap,frame);
            image = ss_convert_stack_to_hue(image);
            image.frameSize = frame * 4;
            bitmap = image.drawBitmapRep(nan);
            imshow(bitmap(:,:,1:3));
        end
    end
    
    delegate.logLoopCount(loop);
    
    regiRec = cat(3,regiRec,regiList.var);
end

if ishandle(messanger)
    delete(messanger);
end
if ishandle(fig)
    delete(fig);
end

result = logical(regiRec);

beep;
end
