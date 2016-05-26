function im = imupsample(oim,varargin)
  % IMUPSAMPLE  Upsample an input image by a factor of 2
  %
  % im = imupsample(oim)
  % im = imupsample(oim,'ParameterName',ParameterValue,...);
  %
  % Inputs:
  %   oim  h by w by c input image
  %   Optional:
  %     'Iters'  followed by number of iterations to 2x upsample {1}
  %     'K'  followed by number of nearest neighbors to interpolate, 1 is
  %       strictly nearest neighbor, >1 uses Shepard weighting {4}
  %     'Radius'  followed by >=1 window radius for matching, {1}
  %     'Sigma'  followed by standard deviation of Gaussian weighting filter,
  %       ~0 is sharp and noisy, >1 tends to be overly smooth {0.8}
  %     'PerChannel'  followed by whether to upsample each of the c channels
  %       independently, this will be faster and sometimes even better {false}
  %     'RGB2NTSC'  followed by whether to convert input RGB image to NTSC
  %       before processing, seems to be a bit better {false}
  % Outputs:
  %   im   2*h by 2*w by c output image
  %

  im=oim;
  % window radius: >=1, but most of the time 1 is good if k>>1
  w = 1;
  % Number of neighbors to interpolate: 1 is nearest neighbor
  k = 4;
  % Controls smoothness: 0<--- sharp and noisy, smooth --->2,inf
  sigma = 0.8;
  % number of 2x upsamples
  iters = 1;
  % important that this defaults to false
  per_channel = false;
  use_ntsc = false;
  
  rec_params = {'Iters','K','Sigma','Radius','RGB2NTSC'};
  rec_vars = {'iters','k','sigma','w','use_ntsc'};
  init_params = {'PerChannel'};
  init_vars = {'per_channel'};
  params = {rec_params{:},init_params{:}};
  vars = {rec_vars{:},init_vars{:}};
  % Map of parameter names to variable names
  params_to_variables = containers.Map(params,vars);
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  % number of channels
  nc = size(im,3);

  if per_channel && nc>1
    im = zeros(size(im,1)*2^iters,size(im,2)*2^iters,nc);
    % Trick to pass on parameters
    rec_vals = cellfun( ...
      @(name)evalin('caller',name),rec_vars,'UniformOutput',false);
    rec_params_vals = reshape({rec_params{:};rec_vals{:}},1,numel(rec_params)*2);
    if use_ntsc
      oim = rgb2ntsc(oim);
    end
    for c = 1:size(oim,3)
      im(:,:,c) = imupsample(oim(:,:,c),rec_params_vals{:});
    end
    if use_ntsc
      im = ntsc2rgb(im);
    end
    return;
  end

  if use_ntsc && nc>1
    im = rgb2ntsc(im);
  end

  for iter = 1:iters
  
    % for every pixel create a 9-long vector of local window
    W = inf([size(im,1) size(im,2) nc (2*w+1)^2]);
    % for every pixel-corner create a 9-long vector of downsampled 9x9 window
    D = inf([size(im,1)-1,size(im,2)-1,nc,(2*w+1)^2]);
    
    % weights
    G = fspecial('gaussian',2*w+1,sigma);
    for x = -w:w
      for y = -w:w
        c = ((w+x)+(2*w+1)*(w+y))+1;
        W(1:end,1:end,:,c) = G(c)*circshift(im(1:end,1:end,:),[x y]);
        D(1:end,1:end,:,c) = G(c)*0.25*( ...
          circshift(im(1:end-1,1:end-1,:),2*[x y]) + ...
          circshift(im(2:end  ,1:end-1,:),2*[x y]) + ...
          circshift(im(2:end  ,2:end  ,:),2*[x y]) + ...
          circshift(im(1:end-1,2:end  ,:),2*[x y]));
      end
    end
    
    % treat every window as a point in R⁹
    W = reshape(W,[],nc*(2*w+1)^2);
    D = reshape(D,[],nc*(2*w+1)^2);
    % find closest match d in D for each w in W
    [I,S] = knnsearch(D,W,'K',k);
    % Shepard inverse distance weights
    S = S.^-2;
    % normalize S weights
    S = bsxfun(@rdivide,S,sum(S,2));
    % If there are perfect matches then there'll be NaNs
    S(isnan(S)) = 1;
    % In case there are duplicates, normalize again
    S = bsxfun(@rdivide,S,sum(S,2));
    
    U = zeros(size(im,1)*2,size(im,2)*2,nc);
    
    for x = [0 1]
      for y = [0 1]
        % Turn index into D into subscript into offset image
        [Y,X] = ind2sub([size(im,1)-1 size(im,2)-1],I);
        Y = min(Y+y,size(im,1));
        X = min(X+x,size(im,2));
        for c = 1:nc
          U( 1+y:2:end , 1+x:2:end,c) = reshape( ...
            sum(im(sub2ind(size(im),Y,X,c*ones(size(X)))).*S,2),...
            size(U,1)/2,size(U,2)/2);
        end
      end
    end
    im = U;
  end

  if use_ntsc && nc>1
    im = ntsc2rgb(im);
  end
  
end
