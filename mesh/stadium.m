function [V,E,F] = stadium(n,a,r)
  % [V,E] = stadium(n,a,r)
  % [V,E,F] = stadium(n,a,r)
  %
  % Construct a polygon approximation of a stadium shape (cigar, pill, capsule)
  % in 2D.
  %
  % Inputs:
  %   n  number of vertices along each semicircle
  %   a  length of straight sides
  %   r  radius
  % Outputs
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into rows of V
  %   F  #F by 3 list of triangle indices into rows of V
  %
  % See also: capsule

  if ~exist('a','var')
    a = 4;
  end
  if ~exist('1','var')
    r = 1;
  end
  assert(a>=0);
  assert(r>=0);

  th = linspace(0,pi,n)';
  h = r*pi/(n-1);
  th = th(2:end-1);

  if a==0
    min_s = 1;
  else
    min_s = 2;
  end
  s = linspace(0.5*a,-0.5*a,max(ceil(a/h),min_s))';
  V = [ ...
    [r*sin(th)+0.5*a -r*cos(th)]; ...
    s ones(numel(s),1); ...
    [-r*sin(th)-0.5*a r*cos(th)]; ...
    flipud(s) -ones(numel(s),1)];
  E = [1:size(V,1);2:size(V,1) 1]';
  if nargout>2
    [V,F] = triangulate(V,E,'Flags',sprintf('-cq32a%0.17f',h^2));
  end
end
