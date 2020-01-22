
function result = ss_color_blue_red(stack)

if isnumeric(stack)
    stack = SYImage(SYData(stack));
else
    stack = stack.copy;
end
if stack.graphicsContext.colorSpace ~= ...
        SYGraphicsContext.ColorSpaceGrayscale
    stack.splitChannels;
end

load('blue-red.lut.mat','lut');
lut(:,4) = 255;
stack.setLut(lut);
stack.range = [-1,1];
stack.graphicsContext.bitsPerComponent = 8;
stack.graphicsContext.colorSpace = SYGraphicsContext.ColorSpaceIndexed;

image = SYImage;
array = stack.bitmapImageArray(true);
for i = 1:array.count
    image.addRepresentation(array.objectAtIndex(i));
end

result = image;
end
