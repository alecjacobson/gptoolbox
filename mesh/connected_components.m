function [C,CF] = connected_components(F)
  % CONNECTED_COMPONENTS Determine the connected components of a mesh described
  % by the simplex list F. Components are determined with respect to the edges of
  % the mesh. That is, a single component may contain non-manifold edges and
  % vertices.
  %
  % C = connected_components(F)
  %
  % Inputs:
  %   F  #F by simplex-size list of simplices
  % Outputs:
  %   C  #V list of ids for each CC 
  %   CF  #F list of ids for each CC
  % 
  % Examples:
  %  trisurf(F,V(:,1),V(:,2),V(:,3), ...
  %    connected_components([F;repmat(size(V,1),1,3)]));

  % build adjacency list
  A = adjacency_matrix(F);
  [~,C] = conncomp(A);
  if nargout > 1 
      CF = C(F(:,1));
  end

end
