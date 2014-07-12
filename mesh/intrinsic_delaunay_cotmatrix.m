function [L,F] = intrinsic_delaunay_cotmatrix(V,F,varargin)
  % INTRINSIC_DELAUNAY_COTMATRIX Construct the intrinsic Delaunay Laplacian,
  % given a mesh (V,F) or (l,F) according to "A discrete Laplace-Beltrami
  % operator for simplicial surfaces" [Bobenko & Springborn 2005] following "An
  % Algorithm for the Construction of Intrinsic Delaunay Triangulations with
  % Applications to Digital Geometry Processing" [Fisher et al. 2006]
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'SideLengths' followed by a #F by 3 list of triangle side lengths
  %       corresponding to edges 23 31 12.
  % Outputs:
  %   L  #V by #V sparse Laplacian matrix
  %   F  #F by 3 list of new intrinsic "triangle" indices, representing
  %     **combinatorially** flipped edges. This is not the "overlay graph" of
  %     [Fisher et al.].
  %

  % default values
  l = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'SideLengths'}, {'l'});
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
  if isempty(l)
    l = edge_lengths(V,F);
  end

  while true
    D = is_delaunay([],F,'SideLengths',l);
    if all(D(:))
      break;
    end
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    NDE = allE(~D,:);
    NDuE = unique(sort(NDE,2),'rows');
    [F,~,l] = flip_edges(F,NDuE,'AllowNonManifold',true,'SideLengths',l,'V',V);
  end

  L = cotmatrix_intrinsic(l,F,size(V,1));

end
