function [V,E] = rounded_rectangle(r,n,w,h,mask)
  % Inputs:
  %   r  radius of corners {0.1}
  %   n  number of samples per corner {16}
  %   w  width {1}
  %   h  height {1}
  %   mask  2by2 mask for corners {ones(2)}
  % Outputs:
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into V
  %

  if nargin<1
    r = 0.1;
  end
  if nargin<2
    n = 16;
  end
  if nargin<3
    w = 1;
  end
  if nargin<4
    h = 1;
  end
  if nargin<5
    mask = ones(2);
  end



  assert(n>1);
  assert(w > 2*r);
  assert(h > 2*r);

  th = linspace(0,pi/2,n)';

  p2c = @(th) [cos(th),sin(th)]*r;

  VV = cell(2,2);

  for ij = [[2;2] [1;2] [1;1] [2;1]]
    [i,j] = deal(ij(1),ij(2));
    if mask(end-j+1,i)
      V{i,j} = [-(2*(i==1)-1)*0.5*(w-2*r) -(2*(j==1)-1)*0.5*(h-2*r)] + p2c(th);
    else
      V{i,j} = [-(2*(i==1)-1)*0.5*(w) -(2*(j==1)-1)*0.5*(h)];
    end
    th = th + pi/2;
  end

  V = [ ...
   V{2,2}
   V{1,2}
   V{1,1}
   V{2,1}
    ];
  E = [1:size(V,1);2:size(V,1) 1]';

end
