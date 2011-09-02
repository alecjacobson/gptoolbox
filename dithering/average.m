% "average" dithering == thresholding on average
av = im > mean(im(:));
% local-average dithering, thresholding on local average, where local average
% is read as blurred original, blur kernel shouldn't be too big. In fact blurry
% is just approximating a blurred segmentated image, which I think would work
% better. If the blur kernel is too big then the result can actually look worse
% than using a constant global average
h = fspecial('gaussian',round(size(im,1)/2),round(size(im,2)/2));
local_average=imfilter(im,h,'replicate');
lav = im > local_average;
% average of global average and local average
lavav = im>(local_average*0.75+0.25*mean(im(:)));
