% Short script to convert a matrix or regions to a stack.
% Written by Satoshi Yamashita.

function result = regiList2image(regiList,data,frame)
% Short script to convert a matrix of regions to an image stack.
% result = regiList2image(regiList,data,frame)
% Argument regiList is a matrix whose rows and columns represent data
% points and regions.
% Argument data is a TSDataMap instance.
% Argument frame is a double[2] representing frame size of the image.
% Return value is an SYImage of grayscale stack whhose slices represent
% regions.

if isa(data,'TSDataMap')
    data = data.dataList;
end

regiN = size(regiList,2);
array = regiList * (1:regiN)';
stack = zeros(frame(1),frame(2),regiN);

for i = 1:length(array)
    if array(i) > 0
        stack(data(i,2),data(i,1),array(i)) = 1;
    end
end

image = SYImage(SYData(stack));
image.splitChannels;

result = image;
end
