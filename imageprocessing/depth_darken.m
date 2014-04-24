function [C,U,Um] = depth_darken(im,D,varargin)
  % DEPTH_DARKEN Enhance an image by darkening according to its corresponding
  % depth image's unsharp mask. "Image Enhancement by Unsharp Masking the Depth
  % Buffer" [Luft et al. 2006].
  %
  % [C,U,Um] = depth_darken(im,D,'ParamterName',ParameterValue, ...)
  %
  % Inputs:
  %   im  h by w by c color image
  %   D  h by w inverse depth image
  %   Optional:
  %     'Lambda' followed by darkening intensity {0.3}
  %     'Sigma' followed by darkening blur radius {0.08*h}
  % Outputs:
  %   C  h by w by c depth-darkened color image 
  %   U  h by w unsharp *mask* of depth image
  %   Um  h by w negative part used for darkening
  %   

  % default values
  lambda = 0.3;
  sigma = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Lambda','Sigma'}, ...
    {'lambda','sigma'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  if isempty(sigma)
    sigma = 0.08*size(D,1);
  end
  % Gaussian blur kernel
  G = fspecial('gaussian',ceil(2*sigma),sigma);
  % Unshapr mask is image minus blurred image
  U = D-imfilter(D,G,'replicate');
  % only keep negative part
  Um = (U<0).*U;
  % Darken according to mask
  C = bsxfun(@plus,im,lambda*Um);

end
