function [T,N,B,I] = per_vertex_frames(V,F,varargin)
  % PER_VERTEX_NORMALS Compute local frames for each vertex of a triangle
  % mesh. This method works by first computing the normal nj at each vertex vi
  % in the mesh (V,F), then finds the neighbor vi of vj such that the
  % undirected line e_ij is most perpendicular to the normal. The cross product
  % of this line defines a binormal and the tangent is recovered
  % by taking the cross product of the binormal and the normal.
  %
  % [T,N,B] = per_vertex_frames(V,F)
  % [T,N,B,I] = per_vertex_frames(V,F,'ParameterName',parameter_value, ...)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'Neighbors' followed by precomputed I (see output below)
  % Outputs:
  %   T  #V by 3 list of tangent vectors
  %   N  #V by 3 list of normal vectors
  %   B  #V by 3 list of binormal vectors
  %   I  #V list of indices of neighbor used to determine normal plane.
  % 
  %

  % default values
  I = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Neighbors'}, ...
    {'I'});
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

  n = size(V,1);
  % normals
  N = per_vertex_normals(V,F);
  if isempty(I)
    % all edges
    E = [F(:,2:3);F(:,[3 1]);F(:,1:2)];
    E = [E;fliplr(E)];
    % All dot products, so that D(i,j) = 1 + n_i ⋅ e_ij
    D = sparse( ...
      E(:,1), ...
      E(:,2), ...
      ... % Important to take max in case e_ij = -n_i
      max( ...
        sum(normalizerow(V(E(:,1),:)-V(E(:,2),:)).*N(E(:,1),:),2), ...
        sum(normalizerow(V(E(:,2),:)-V(E(:,1),:)).*N(E(:,1),:),2)) ...
      +1, ...
      n,n);
    [~,I] = minnz(D);
  end
  % edge-vector that looks most like a tangent vector
  Te = normalizerow(V(I,:)-V);
  % Bi-normal orthogonal to normal-plane cutting most tangent edge-vector 
  B = normalizerow(cross(Te,N,2));
  % Finally recompute a tangent vector orthogonal to normal and binormal
  T = cross(B,N,2);

end
