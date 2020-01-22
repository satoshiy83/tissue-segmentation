% Short script drawing a bitmap with silhouette value.
% Written by Satoshi Yamashita.

function result = ss_draw_silhouette_map(partition,data,hint)
% Short script to draw a bitmap with silhouette value.
% result = ss_draw_silhouette_map(partition,data,hint)
% Argument partition is a matrix whose rows and columns represent data
% points and regions.
% Argument data is a TSDataMap instance.
% Argument hint is an SYDictionary instance holding parameters for
% silhouette analysis. For parameters, see TSMeter.
% Return value is a grayscale SYImage instance.

analyser = TSRegionAnalyser(data);
meter = TSMeter(data,hint);
analyser.meter = meter;
partData = SYData(partition);
analyser.regiList = partData;
meter.regiList = partData;
% analyser.regiN = size(partition,2);

silhList = analyser.silhouette;
silhList = silhList';
silhList = silhList(partition');
frame = data.frameSize;
bitmap = nan(frame,'single');

for i = 1:length(silhList)
    bitmap(data.dataList(i,2),data.dataList(i,1)) = silhList(i);
end

% image = SYImage(bitmap);
context = SYGraphicsContext(SYData(bitmap),frame(2),frame(1),32, ...
    SYGraphicsContext.CompositeModeOver, ...
    SYGraphicsContext.ColorSpaceGrayscale, nan, [-1,1]);
image = SYImage(context);

result = image;
end
