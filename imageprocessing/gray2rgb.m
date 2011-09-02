function rgb = gray2rgb(im,map,n)
  % GRAY2RGB
  %
  % rgb = gray2rgb(im,map,n)
  %
  % Use a colormap, map, to convert a grayscale/intensity image with n colors
  % to pseudocolor rgb image. The n colors in the original image are mapped
  % evenly to colors in the given colormap, not according to value.
  %
  % rgb = gray2rgb(im,map)
  %
  % Same as above but compute number of colors in im
  %
  % rgb = gray2rgb(im)
  %
  % Use default colormap and compute number of colors in im
  %
  % Inputs:
  %  im  grayscale image, 2d array (h,w)
  %  map  colormap
  %  n  (optional) number of colors in im, default is to compute unique colors,
  %    this is usually preferred but it's time consuming so if the value is
  %    already known then it can be given.
  % Outputs:
  %  rgb  color Red, Green, Blue image, 3d array (h,w,3)
  %
  % See also: rgb2gray, hsv2rgb, rgb2hsv, imshow
  %

  % number of unique colors
  if(~exist('n','var'))
    n = size(unique(im),1);
  end

  if(~exist('map','var'))
    if(~exist('gray2ind'))
      error( ...
        'Image processing toolbox seems not installed, try again without map')
    end
    if(isempty(get(0,'CurrentFigure')))
      map = jet(n);
    else
      map = colormap;
    end
    %map = jet(n);
    m = size(map,1);
  else
    m = size(map,1);
  end

  if(exist('gray2ind')&&n<=65536)
    
    in = gray2ind(im,n);

    if(strcmp(class(in),'uint8') && (n/(n/m)) > 256)
      in = uint16(in);
    end

    if( m <= 256 )
      rgb = label2rgb(idivide(in,n/m,'fix')+1,map);
    else
      rgb = ind2rgb(idivide(in,n/m,'fix'),map);
    end

    % if it's not going to change the result it's faster to convert to and
    % from uint8 and use label2rgb
    if(strcmp(class(im),'double'))
      rgb = im2double(rgb);
    elseif(strcmp(class(im),'uint16'))
      rgb = im2uint16(rgb);
    end
  else 
    % slope
    m = 1/(1/4);
    % precomputation
    mim = m*im;
    nmim = -mim;
    x1 = im<1/8;
    x2 = im>=1/8 & im<3/8;
    x3 = im>=3/8 & im<5/8;
    x4 = im>=5/8 & im<7/8;
    x5 = im>=7/8;
    rgb = cat( 3, ...
      x3.*(mim - (3/8)*m)   + x4.*1 + x5.*(nmim + (9/8)*m), ...
      x2.*(mim - (1/8)*m)   + x3.*1 + x4.*(nmim + (7/8)*m), ...
      x1.*(mim + 0.5)       + x2.*1 + x3.*(nmim + (5/8)*m));
  end
end
