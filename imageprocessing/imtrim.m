function [C,top, bottom, left, right] = imtrim(im,varargin)
  % IMTRIM auto-crop an image like Photoshop's Edit>Trim feature, as of yet only
  % grayscale imagse are supported
  %
  % C = imtrim(im)
  % [C,top,bottom,left,right] = imtrim(im,'ParamName',ParamValue, ...)
  %   
  % Inputs:
  %   im  h by w by nc image
  %   Optional:
  %     'Location'  followed by {'NorthWest'} or 'SouthEast' specifying where
  %       to find trim color
  %     'Dims' followed by 2 long logic whether to use each dimension: [1 0]
  %       means only crop in h and not in w
  %     'Threshold' followed by distance threshold {0.1}
  %     'Map'  followed by colormap {[]}
  % Outputs:
  %   C  sh by wh cropped image
  %   top  index of top pixel in im
  %   bottom  index of bottom pixel in im
  %   left  index of left pixel in im
  %   righ  index of right pixel in im
  %
  %   Copyright Alec Jacobson, 2010
  %

  location = 'NorthWest';
  threshold = [];
  dims = [1 1];
  map = [];

  % parse optional input parameters
  v = 1;
  while v < numel(varargin)
    switch varargin{v}
    case 'Map'
      assert(v+1<=numel(varargin));
      v = v+1;
      map = varargin{v};
    case 'Location'
      assert(v+1<=numel(varargin));
      v = v+1;
      location = varargin{v};
    case 'Threshold'
      assert(v+1<=numel(varargin));
      v = v+1;
      threshold = varargin{v};
    case 'Dims'
      assert(v+1<=numel(varargin));
      v = v+1;
      dims = varargin{v};
    otherwise
      error('Unsupported parameter: %s',varargin{v});
    end
    v = v+1;
  end

  if isempty(threshold)
    % default threshold parameter, works equivalently with Photoshop's
    % hardcoded parameter
    threshold = 0.1;
  end

  im_in = im;

  if ~isempty(map)
    nf = size(im,4);
    im = cell2mat(permute(arrayfun(@(C) ind2rgb(im(:,:,:,C),map(:,:,min(end,C))),1:size(im,4),'UniformOutput',false),[1 4 3 2]));
  end

  if isempty(map) && ~isfloat(im)
    warning('converting input to double');
    im = im2double(im);
  end

  % gather corner value to which the image is compared
  if(strcmp(location,'NorthWest'))
    corner_value = im(1,1,:,1);
  elseif(strcmp(location,'SouthEast'))
    corner_value = im(1,1,:,1);
  else
    error([location ' is not a valid location']);
  end

  % get difference of image with corner value
  %difference = abs(im - corner_value)>0.1;
  % should work for any number of channels
  %difference = sqrt(sum((im - corner_value).^2,3)) > ...
  %  sqrt(threshold^2*size(im,3)); 
  difference = sqrt(sum(bsxfun(@minus,im,corner_value).^2,3)) > ...
    sqrt(threshold^2*size(im,3)); 
  difference = max(difference,[],4);

  if dims(2)
    [left_i,left] = ind2sub(size(difference),find(difference,1));
    [right_i,right] = ind2sub(size(difference),find(difference,1,'last'));
  else
    left = 1;
    right = size(im,2);
  end
  if dims(1)
    [top_j,top] = ind2sub(size(difference'),find(difference',1));
    [bottom_j,bottom] = ind2sub(size(difference'),find(difference',1,'last'));
  else
    top = 1;
    bottom = size(im,1);
  end
  C = im_in(top:bottom,left:right,:,:);
end
