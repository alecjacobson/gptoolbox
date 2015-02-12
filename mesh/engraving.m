function [V,F] = engraving(im,w,t,s,l)
  % ENGRAVING Create an engraving of an image
  %
  % [V,F] = engraving(im,w,t,s)
  %
  % Inputs:
  %   im  width by height grayscale image
  %   w   desired width of engraving model (in mm)
  %   t  desired thickness of engraving model (in mm)
  %   s  desired span of thickness devoted to levels (in mm)
  %   l  number of levels
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %
  % Example:
  %  % load image
  %  im = imresize(rgb2gray(im2double(imread('hans-hass.jpg'))),0.5);
  %  % pad image by 10% of width
  %  im = padarray(im,ceil(0.1*repmat(size(im,2),1,2)));
  %  % engrave: 50mm wide, 5mm thick, 1mm devoted to 4 layers
  %  [V,F] = engraving(im,50,5,1,4);
  %

  assert(isfloat(im));
  [V,F] = create_regular_grid(size(im,2),size(im,1),0,0);
  V(:,1) = V(:,1)*w;
  V(:,2) = (1-V(:,2))*size(im,1)/size(im,2)*w;
  [V,F] = extrude(V,F);
  V(:,3) = V(:,3)*t;
  V(1:numel(im),3) = V(1:numel(im),3)-s*round(im(:)*(l-1))/(l-1);

end
