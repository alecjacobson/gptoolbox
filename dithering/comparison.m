im = im2double(imread('david.png'));
% threshold
th = round(im);
% random dither
ra = im > rand(size(im));
% random dither, threshold random choice
rr = round(rand(size(im)));
rara = rr.*round(im)+(1-rr).*(im > rand(size(im)));
% high contrast, random dither
sp = -(im).*(im).*(im-1.5)*2>rand(size(im));
% random dither, threshold, random choice except keep edges
e = edge(im,0.05);
blur_width=5;
h = fspecial('gaussian',round(blur_width),round(blur_width));
blurred_e=ceil(imfilter(e+0.0,h,'replicate'));
rre = (1-((1-rr).*(1-blurred_e)));
rarae = rre.*round(im)+(1-rre).*(im > rand(size(im)));
imwrite(im,'dithering-original.png');
imwrite(th,'dithering-thresh.png');
imwrite(ra,'dithering-random.png');
imwrite(rara,'dithering-random-thresh.png');
imwrite(sp,'dithering-random-contrast.png');
imwrite(rarae,'dithering-random-thresh-edges.png');
