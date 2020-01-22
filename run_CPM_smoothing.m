% Smoothing of segments boundary by Cellular Potts model.
% Written by Satoshi Yamashita.
% 
% This function run a smoothing of segments boundary by Cellular Potts
% model with Hamiltonian of area, interface, and silhouette.

function result = run_CPM_smoothing(partition,data,hint,delegate)
% Function implementing cellular Potts model simulation with silhouette.
% result = run_CPM_smoothing(partition,data,hint,delegate)
% Argument partition is a matrix whose rows and columns represent data
% points and regions.
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters.
% Parameters:
%   cpm_iterN: (double) maximum number to try updates.
%   cpm_uplaN: (double) number to update pixel label.
%   cpm_coefficients: (double[1,3]) array of coefficients for area
%   constraint, surface tension, and homogeneity of region.
%   cpm_temperature: (double) tempareture of system.
% Argument delegate is a CPMDelegate instance.
% Return value is a logical matrix whose rows and columns represent data
% points and clusters.

%% Initialization.
iterN = 20000;
if ~isnan(hint.objectForKey('cpm_iterN'))
    iterN = hint.objectForKey('cpm_iterN');
end
uplaN = 100;
if ~isnan(hint.objectForKey('cpm_uplaN'))
    uplaN = hint.objectForKey('cpm_uplaN');
end

regiData = SYData(partition);
meter = TSMeter(data,hint);
meter.regiList = regiData;
analyser = TSRegionAnalyser(data);
analyser.regiList = regiData;
analyser.meter = meter;

dict = hint.copy;
areas = analyser.size_regions;
dict.setObjectForKey('areas',areas);

coef = hint.objectForKey('cpm_coefficients');
Hamilements = {@lf_area, @lf_surface, @lf_silhouette};
if any(isnan(coef)) || (length(coef) ~= length(Hamilements))
    disp('Hint must contain ''cpm_coefficents'' of correct length!');
    result = nan;
    return;
end
H = 0;
for k = 1:length(Hamilements)
    fh_Hamilement = Hamilements{k};
    H = H + fh_Hamilement(analyser,dict) * coef(k);
end

T = 1;
if ~isnan(hint.objectForKey('cpm_temperature'))
    T = hint.objectForKey('cpm_temperature');
end

%% Main loop: change of label for a point is tested in each step of t.
m_dH = 0; c_ndH = 0; uplaC = 0; count = 1;
messanger = dialog('Position',[300 300 200 150],'WindowStyle','normal', ...
    'Name','Cellular Potts model');
msg1 = uicontrol(messanger,'Style','text','Position',[5 35 190 40],'FontSize',14);
msg2 = uicontrol(messanger,'Style','text','Position',[5 85 190 40],'FontSize',14);
button = uicontrol(messanger,'Style','togglebutton','Position',[55 5 90 30], ...
    'FontSize',14,'String','Exit');
for t = 1:iterN
    % User interface.
    count = count + 1;
    if count > 10
        count = 1;
        if ishandle(messanger)
           msg1.String = ['Iterating ',num2str(t),' / ',num2str(iterN)];
           drawnow;
        end
    end
    if ishandle(messanger)
        if button.Value
            break;
        end
    end
    
    % Select a boundary point randomly.
    rim = analyser.rim_points;
    rim = find(rim);
    if isempty(rim)
        break;
    end
    i = ceil(rand() * length(rim));
    i = rim(i);
    
    % Select a label from neighbors of the point.
    neighbors = data.neighborsOfPoint(i);
    neighbors = neighbors(randperm(length(neighbors)));
    j = neighbors(1);
    if j < 1 || any(regiData.var(i,:) & regiData.var(j,:))
        continue;
    end
    
    % Check connectedness.
    m = zeros(3,3,'logical');
    x_0 = data.dataList(i,1) - 2;
    y_0 = data.dataList(i,2) - 2;
    for k = neighbors
        x = data.dataList(k,1) - x_0;
        y = data.dataList(k,2) - y_0;
        m(y,x) = any(regiData.var(i,:) & regiData.var(k,:));
    end
    array = m([1,2,3,6,9,8,7,4]);
    brray = m([2,3,6,9,8,7,4,1]);
    if sum(array ~= brray) > 2
        continue;
    end
    
    % Get a change of energy for a change of label.
    oldLabel = regiData.var(i,:);
    regiData.var(i,:) = regiData.var(j,:);
    I = 0;
    for k = 1:length(Hamilements)
        fh_Hamilement = Hamilements{k};
        I = I + fh_Hamilement(analyser,dict) * coef(k);
    end
    dH = I - H;
    
    % Adapt the change of label on the probability.
    p = exp(-dH / T);
    if p > rand()
        m_dH = m_dH + dH;
        if dH < 0
            c_ndH = c_ndH + 1; 
        end
        if ishandle(messanger)
            msg2.String = ['Count ',num2str(uplaC),' / ',num2str(uplaN), ...
                '; dH: ',num2str(dH,'%.4g'),'; p: ',num2str(p,'%.4g')];
            drawnow;
        end
        
        H = I;
        
        uplaC = uplaC + 1;
        if uplaC >= uplaN
            break;
        end
    else
        regiData.var(i,:) = oldLabel;
    end
end
m_dH = m_dH / uplaN;
disp(['Mean energy change: ',num2str(m_dH),'; Count of enegy decreasing steps :',num2str(c_ndH),'/',num2str(uplaC)]);
delegate.record = SYDictionary({'m_dH',m_dH; 'c_ndH',c_ndH; 'uplaC',uplaC});

if ishandle(messanger)
    delete(messanger);
end

result = regiData.var(:,any(regiData.var,1));

end

function result = lf_area(analyser,hint)
areas = hint.objectForKey('areas');
breas = analyser.size_regions;

result = sum((areas - breas) .^ 2);
end
function result = lf_surface(analyser,~)
interfaces = analyser.count_connections;

result = sum(interfaces);
end
function result = lf_silhouette(analyser,~)
sil = analyser.silhouette;

result = size(sil,1) - sum(sum(sil));
end
function result = lf_chechConnectedness(analyser,hint)

end
