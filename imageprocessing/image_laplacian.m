function L = image_laplacian(varargin)
  % IMAGE_LAPLACIAN Compute the image laplacian for a given image in the manner
  % of "Colorization using Optimization" by [Levin et al. 2004]. This is a sort
  % of amalgamation of the literal description in the paper but the knowledge
  % that a Laplacian (rather than a bi-Laplacian) is being used.
  % 
  % L = IMAGE_LAPLACIAN(im)
  % L = IMAGE_LAPLACIAN(im,'ParameterName',ParameterValue)
  %
  % Inputs:
  %   im  h by w by (3|1) image (expects double)
  %   Optional:
  %     'Omega' followed by an omega value
  % Outputs:
  %   L w*h by w*h sparse laplacian matrix
  %
  % See also: cotmatrix, levin, lischinski
  %

  im = varargin{1};

  % width and height
  h = size(im,1);
  w = size(im,2);
  %   image size
  %  w     h    c  harmonic biharmonic
  %  185   138  1  0.05     0.1
  %  185   138  3           0.1
  %  93    69   3           0.1
  omega = 0.1;
  % number of vertices/pixels
  n = h*w;

  if ~strcmp(class(im),'double')
    warning('Casting input to double: `im = im2double(im)`');
    im = im2double(im);
  end

  ii = 2;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Omega'
      ii = ii + 1;
      assert(ii<=nargin)
      omega = varargin{ii};
    otherwise
      error('Unsupported parameter');
    end
    ii = ii+1;
  end

  % runs across height fastest
  I = (1:n)';
  I = reshape(I,h,w);
  Ileft = I(1:h,1:(w-1));
  Iright = I(1:h,2:w);
  Ibottom = I(1:(h-1),1:w);
  Itop= I(2:h,1:w);
  % determine neighborhood via FD stencil
  E = [ ...
    Ileft(:) Iright(:);...
    Iright(:) Ileft(:);...
    Ibottom(:) Itop(:);...
    Itop(:) Ibottom(:);...
    ];
  assert(size(E,1) == ((h-1)*w + (w-1)*h)*2);
  r = E(:,1);
  s = E(:,2);
  imlong = reshape(im,size(im,1)*size(im,2),size(im,3));

  wrs = exp(-sum((imlong(r,:)-imlong(s,:)).^2,2)./(2*omega.^2));

  L = sparse(r,s,wrs,n,n);
  % We've made L(r,s) = wrs, now we need L(r,r) = ∑L(r,s) and L(r,s) = -wrs
  L = diag(sum(L,2))-L;
end
