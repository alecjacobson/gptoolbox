function [D] = div(V,F)
  % DIV Compute the numerical divergence at each vertex of a mesh of a vector
  % field defined at each face of a mesh at every face of a triangle mesh.
  %
  % D = div(V,F)
  %
  % Inputs:
  %   V  #vertices by 3 list of mesh vertex positions
  %   F  #faces by 3 list of mesh face indices
  % Outputs:
  %   D  #vertices by #faces*3 divergence operator
  %

  %% append with 0s for convenience
  %if size(V,2) == 2
  %  V = [V zeros(size(V,1),1)];
  %  if size(X,2) == 2
  %    X = [X zeros(size(X,1),1)];
  %  end
  %end

  %% number of domain vertices 
  %n = size(V,1);

  %assert(size(V,2) == size(X,2));
  %assert(size(X,1) == size(F,1));

  %% renaming indices of vertices of triangles for convenience
  %i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
  %% #F x 3 matrices of triangle edge vectors, named after opposite vertices
  %v23 = V(i3,:) - V(i2,:);  v31 = V(i1,:) - V(i3,:); v12 = V(i2,:) - V(i1,:);

  %% The integrated divergence associated with vertex 𝑖 can be written as
  %% ∇ · 𝑋 = 1 ∑︁ cot 𝜃1 (𝑒1 · 𝑋𝑗 ) + cot 𝜃2 (𝑒2 · 𝑋𝑗 )

  %C = cotangent(V,F);
  %DX = sparse( ...
  %  [i1;i2;i3], ...
  %  1, ...
  %  [ ...
  %    C(:,2) .* sum(-v31.*X,2) + C(:,3).*sum(v12.*X,2) ;...
  %    C(:,3) .* sum(-v12.*X,2) + C(:,1).*sum(v23.*X,2) ;...
  %    C(:,1) .* sum(-v23.*X,2) + C(:,2).*sum(v31.*X,2) ;...
  %  ], ...
  %  n,1);

  %DX = full(0.5*DX);

  % This seems to be identical to:
  G = grad(V,F); 
  switch size(F,2)
  case 3
    dblvol = doublearea(V,F);
  case 4
    dblvol = 2*volume(V,F);
  end
  TA = repdiag(diag(sparse(dblvol)),size(V,2));
  D = -0.25*G'*TA;
  %DX = D*X(:);
end
