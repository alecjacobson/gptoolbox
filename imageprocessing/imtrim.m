function [out1,out2,out3,out4,out5] = imtrim(im,location)
  %IMTRIM auto-crop an image like Photoshop's Edit>Trim feature, as of yet only
  %grayscale imagse are supported
  %
  %   cropped = IMTRIM(IM) crop image based on top left corner
  %   
  %   [cropped,t,b,l,r] = IMTRIM(IM) return cropped image and indices used to
  %   crop the image. So cropped = im(t:b,l:r);
  %
  %   [t,b,l,r] = IMTRIM(IM) return only indices used to crop the image. So 
  %   cropped = im(t:b,l:r);
  %
  %   [...] = IMTRIM(IM,location) same as above but location may specify
  %   top-left corner ('NorthWest') or bottom-right corner ('SouthEast') to be
  %   the picel used in determining the auto-crop
  %
  %   Copyright Alec Jacobson, 2010
  %

  if(~exist('location'))
    location = 'NorthWest';
  end

  % gather corner value to which the image is compared
  if(strcmp(location,'NorthWest'))
    corner_value = im(1,1);
  elseif(strcmp(location,'SouthEast'))
    corner_value = im(1,1);
  else
    error([location ' is not a valid location']);
  end

  % hard-coded threshold parameter, works equivalently with Photoshop's
  % hardcoded parameter
  threshold = 0.1;

  % get difference of image with corner value
  %difference = abs(im - corner_value)>0.1;
  % should work for any number of channels
  difference = sqrt(sum((im - corner_value).^2,3)) > ...
    sqrt(threshold^2*size(im,3)); 
  [left_i,left_j] = ind2sub(size(difference),find(difference,1));
  [right_i,right_j] = ind2sub(size(difference),find(difference,1,'last'));
  [top_j,top_i] = ind2sub(size(difference'),find(difference',1));
  [bottom_j,bottom_i] = ind2sub(size(difference'),find(difference',1,'last'));
  if(nargout == 1)
    out1 = im(top_i:bottom_i,left_j:right_j);
  elseif(nargout == 5)
    out1 = im(top_i:bottom_i,left_j:right_j);
    out2 = top_i;
    out3 = bottom_i;
    out4 = left_j;
    out5 = right_j;
  else
    out1 = top_i;
    out2 = bottom_i;
    out3 = left_j;
    out4 = right_j;
  end
end
