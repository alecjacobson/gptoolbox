function V = vhs(im,varargin)
  % VHS Apply a VHS filter to an image
  %
  % V = vhs(im)
  % V = vhs(im,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   im  h by w by c image
  %   Optional:
  %     'VerticalLoop' followed by whether to loop the bent strip vertically
  %     over time and output a sequence of images.
  % Output:
  %   V  h by w by c by f image of (f long sequence of images)
  %

  % must use double
  im = im2double(im);

  looping = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'VerticalLoop'},{'looping'});
  v = 1;
  iter = 1;
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

  if looping 
    strip_top = 1;
  else
    strip_top = ceil(0.4*size(im,1));
  end

  first = 1;
  while true

    % http://mikeyjam.buzznet.com/user/journal/12237761/tutorial-getting-vhs-tv-effect/
    A = im;
    B = im;
    C = im;
    D = im;
    % exclusion blend like photoshop
    % http://www.deepskycolors.com/archive/2010/04/21/formulas-for-Photoshop-blending-modes.html
    ex = @(T,B) 0.5 - 2.*bsxfun(@times,T-0.5,B-0.5);
    % kill color channels
    A(:,:,1) = 0;
    B(:,:,2) = 0;
    C(:,:,3) = 0;
    % Shift color layers
    nudge = @(f) ceil(f*rand(1)*size(im,2));
    A = A(:,mod(nudge(0.02)+(1:end)-1,end)+1,:);
    B = B(:,mod(nudge(0.02)+(1:end)-1,end)+1,:);
    C = C(:,mod(nudge(0.02)+(1:end)-1,end)+1,:);
    A = A(mod(nudge(0.005)+(1:end)-1,end)+1,:,:);
    B = B(mod(nudge(0.005)+(1:end)-1,end)+1,:,:);
    C = C(mod(nudge(0.005)+(1:end)-1,end)+1,:,:);
    % exclusion blend colored layers and alpha blend with original
    
    F = D+0.3*(ex(ex(C,B),A)-D);
    N = rand(size(im));

    % inverse mapping function
    bend_w = (2*rand(1)-1)*5;
    bend = @(x,u) [mod(x(:,1)+(1-x(:,2)/max(x(:,2))).^2*bend_w,max(x(:,1))) x(:,2)];
    % maketform arguments
    ndims_in = 2;
    ndims_out = 2;
    tdata = [];
    tform = maketform('custom', 2,2, [], bend, tdata);

    % Bend strip
    strip_h = ceil(1/4*size(im,1));
    strip = mod(strip_top+(1:strip_h)-1,size(im,1))+1;
    F(strip,:,:) = imtransform(F(strip,:,:), tform);

    % overlay gray line
    ol= @(T,B) (T>0.5).*(1-(1-2.*(T-0.5)).*(1-B))+(T<=0.5).*((2.*T).*B);
    G = repmat(0.75,[numel(-1:1) size(F,2) size(F,3)]);
    F(mod(strip_top+(-1:1)-1,size(F,1))+1,:,:) = ...
      ol( F(mod(strip_top+(-1:1)-1,size(F,1))+1,:,:),G);


    % overlay horizontal lines
    L = zeros(size(F));
    L(1:4:end,:,:) = 1;
    L(2:4:end,:,:) = 1;
    L = imfilter(L,fspecial('gaussian',[5 5],1.5),'replicate');
    F = ol(F,L);

    % Fade in random color gradient
    m = rand(1,2)*2-1;
    R = bsxfun(@plus,m(1)*(1:size(F,2)),m(2)*(1:size(F,1))');
    R = (R-min(R(:)))./(max(R(:))-min(R(:)));
    R = gray2rgb(R,colormap(jet(255)));
    F = F+0.05*(ol(F,R)-F);

    % soft light
    sl = @(T,B) (B>0.5).*(1-(1-T).*(1-(B-0.5))) + (B<=0.5).*(T.*(B+0.5));
    F = F + 0.15*(sl(F,N)-F);

    % sharpen
    V(:,:,:,iter) = imsharpen(F);

    %imshow(V(:,:,:,iter));
    %drawnow;

    strip_top = strip_top+round(0.057*size(im,1));
    iter = iter+1;
    if ~looping || strip_top > size(im,1)
      break;
    end

  end
end
