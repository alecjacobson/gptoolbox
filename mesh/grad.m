function [G] = grad(V,F)
  % GRAD
  % g = grad(V,F)
  %
  % Compute the numerical gradient operator
  %
  % Inputs:
  %   V  #vertices by 3 list of mesh vertex positions
  %   F  #faces by 3 list of mesh face indices
  % Outputs:
  %   G  #faces*3 by #V Gradient operator
  %
  % Example:
  %   L = cotmatrix(V,F)
  %   G  = grad(V,F);
  %   dblA = doublearea(V,F);
  %   GMG = -G'*repdiag(diag(sparse(dblA)/2),3)*G;
  %

  % append with 0s for convenience
  dim = size(V,2);
  if size(V,2) == 2
    V = [V zeros(size(V,1),1)];
  end

  % Gradient of a scalar function defined on piecewise linear elements (mesh)
  % is constant on each triangle i,j,k:
  % grad(Xijk) = (Xj-Xi) * (Vi - Vk)^R90 / 2A + (Xk-Xi) * (Vj - Vi)^R90 / 2A
  % where Xi is the scalar value at vertex i, Vi is the 3D position of vertex
  % i, and A is the area of triangle (i,j,k). ^R90 represent a rotation of 
  % 90 degrees
  %
  % renaming indices of vertices of triangles for convenience
  i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
  % #F x 3 matrices of triangle edge vectors, named after opposite vertices
  v32 = V(i3,:) - V(i2,:);  v13 = V(i1,:) - V(i3,:); v21 = V(i2,:) - V(i1,:);

  % area of parallelogram is twice area of triangle
  % area of parallelogram is || v1 x v2 || 
  n  = cross(v32,v13,2); 
  
  % This does correct l2 norm of rows, so that it contains #F list of twice
  % triangle areas
  dblA = normrow(n);
  
  % now normalize normals to get unit normals
  u = normalizerow(n);

  % rotate each vector 90 degrees around normal
  eperp21 = bsxfun(@times,normalizerow(cross(u,v21)),normrow(v21)./dblA);
  eperp13 = bsxfun(@times,normalizerow(cross(u,v13)),normrow(v13)./dblA);

  %g =                                              ...
  %  (                                              ...
  %    repmat(X(F(:,2)) - X(F(:,1)),1,3).*eperp13 + ...
  %    repmat(X(F(:,3)) - X(F(:,1)),1,3).*eperp21   ...
  %  );

  G = sparse( ...
    [0*size(F,1)+repmat(1:size(F,1),1,4) ...
    1*size(F,1)+repmat(1:size(F,1),1,4) ...
    2*size(F,1)+repmat(1:size(F,1),1,4)]', ...
    repmat([F(:,2);F(:,1);F(:,3);F(:,1)],3,1), ...
    [eperp13(:,1);-eperp13(:,1);eperp21(:,1);-eperp21(:,1); ...
     eperp13(:,2);-eperp13(:,2);eperp21(:,2);-eperp21(:,2); ...
     eperp13(:,3);-eperp13(:,3);eperp21(:,3);-eperp21(:,3)], ...
    3*size(F,1), size(V,1));

  if dim == 2
    G = G(1:(size(F,1)*dim),:);
  end


  % Should be the same as:
  % g = ... 
  %   bsxfun(@times,X(F(:,1)),cross(u,v32)) + ...
  %   bsxfun(@times,X(F(:,2)),cross(u,v13)) + ...
  %   bsxfun(@times,X(F(:,3)),cross(u,v21));
  % g = bsxfun(@rdivide,g,dblA);
end
