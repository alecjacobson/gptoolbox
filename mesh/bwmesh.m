function [W,F,V,E,H] = bwmesh(A,varargin)
  % BWMESH Construct a mesh from a black and white image.
  %
  % [W,F] = bwmesh(A)
  % [W,F] = bwmesh(png_filename)
  % [W,F,V,E,H] = bwmesh(...,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   A  a h by w black and white image (grayscale double images will be
  %     clamped)
  %     or
  %   png_filename  path to a png file
  %   Optional:
  %     'Tol'  tolerance for Douglas-Peucker algorithm {0}
  %     'TriangleFlags' followed by flags to pass to triangle
  %       {'-q30a[median_sqr_edge_length]'}
  %     'SmoothingIters' followed by smoothing amount {0}
  % Outputs:
  %   W  #W by 2 list of mesh vertices
  %   F  #F by 3 list of triangle indices into W
  %   V  #V by 2 list of boundary polygon vertices
  %   E  #E by 2 list of boundary edge indices into V
  %   H  #H by 2 list of hole indicator point positions
  %
  % Known issues: 
  %   - this does _not_ trace the boundary of pixels but rather
  %     connects the centers of boundary pixels. Therefore it struggles if
  %     there are very thin (1px wide) features/holes. 
  %   - This will correctly carve out simple holes, but will not find islands
  %     in holes: e.g. a bull's eye sign.
  %

  if ischar(A)
    % read alpha channel
    [~,~,A] = imread(A);
  end
  A = im2double(A);

  % default values
  tol = 0;
  triangle_flags = '';
  smoothing_iters = 0;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol','TriangleFlags','SmoothingIters'}, ...
    {'tol','triangle_flags','smoothing_iters'});
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

  % B contains list of boundary loops then hole loops, N number of outer
  % boundaries (as opposed to hole boundaries)
  [B,~,N] = bwboundaries(A>0.5);
  % If you don't have the image processing toolbox you can use this (slower)
  % version:
  %[B,~,N] = gp_bwboundaries(A>0.5);

  V = [];
  E = [];
  H = [];
  for b = 1:numel(B)
    if size(B{b},1)>2
      Vb = B{b};
      Vb = bsxfun(@plus,Vb*[0 -1;1 0],[-0.5,size(A,1)+0.5]);
      if smoothing_iters > 0
        Eb = [1:size(Vb,1);2:size(Vb,1) 1]';
        Vb = curve_smooth(Vb,Eb,'MaxIters',smoothing_iters);
      end
    
      if tol > 0 
        Vb = dpsimplify(Vb([1:end 1],:),tol);
        Vb = Vb(1:end-1,:);
      end
      % don't consider degenerate boundaries
      area = sum(prod(Vb([2:end 1],:)+Vb*[-1 0;0 1],2))/2;
      if size(unique(Vb,'rows'),1)>2 && area>eps
        Eb = [1:size(Vb,1);2:size(Vb,1) 1]';
        E = [E;size(V,1)+Eb];
        V = [V;Vb];
        if b > N
          H = [H;point_inside_polygon(Vb)];
        end
      end
    end
  end
  %fprintf('triangle...\n');
  if isempty(triangle_flags)
    % triangulate the polygon
    % get average squared edge length as a guess at the maximum area constraint
    % for the triangulation
    median_sqr_edge_length = median(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2))/2.0;
    quality = 30;
    triangle_flags = sprintf('-q%da%0.17f',quality,median_sqr_edge_length);
  end

  % Why do I need to do this?
  [V,~,J] = remove_duplicate_vertices(V,0);
  % remap faces
  E = J(E);

  [W,F] = triangle(V,E,H,'Flags',triangle_flags,'Quiet');
end
