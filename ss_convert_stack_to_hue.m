% Short script to convert a stack of images to a colored image.
% Written by Satoshi Yamashita.

function result = ss_convert_stack_to_hue(image)
% Short script to convert a stack to a colored image.
% result = ss_convert_stack_to_hue(image)
% Argument image is an SYImage instance of grayscale.
% Return value is an SYImage instance of RGBA in which slices were given a
% color.

frame = image.frameSize;
bitmap = zeros(frame(1),frame(2),4,'uint8');
rgb = ones(1,4); rgb(4) = 255;

% stack = logical(image.stackImageData);
% count = image.countStack;

array = image.bitmapImageArray(false);
count = array.count;
for h = 1:count
%     color = SYColor([h / count,1,0.5],SYColor.ColorHSB);
%     rgb(1:3) = color.RGBColor * 255;
    color = SYColor(8,SYColor.ColorHSB, ...
        uint8([h / count,1,0.5] * 255),nan,nan);
    rgb(1:3) = color.RGBColor;
    
%     mask = stack(:,:,1,h);
    bitmapRep = array.objectAtIndex(h);
    mask = logical(bitmapRep.bitmap.var);
    
    bitmap(mask(:,:,[1,1,1,1])) = repmat(uint8(rgb),sum(mask(:)),1);
end

result = SYImage(SYData(bitmap));
end
