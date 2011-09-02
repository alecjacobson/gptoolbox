function G = grad(V,F,X)
  % GRAD
  % G = grad(V,F,X)
  %
  % Compute the numerical gradient at every face of a triangle mesh.
  %
  % Inputs:
  %   V  #vertices by 3 list of mesh vertex positions
  %   F  #faces by 3 list of mesh face indices
  %   X  # vertices list of scalar function values
  % Outputs:
  %   G  #faces by 3 list of gradient values
  %

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
  eperp21 = normalizerow(cross(u,v21)) .* repmat(normrow(v21),1,3);
  eperp13 = normalizerow(cross(u,v13)) .* repmat(normrow(v13),1,3);

  G =                                              ...
    (                                              ...
      repmat(X(F(:,2)) - X(F(:,1)),1,3).*eperp13 + ...
      repmat(X(F(:,3)) - X(F(:,1)),1,3).*eperp21   ...
    )                                              ...
    ./repmat(dblA,1,3);
end
