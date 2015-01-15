function L = cotmatrix(V,F)
  % COTMATRIX computes cotangent matrix (laplacian mesh operator), (mass/area
  % terms already cancelled out)
  %
  % L = cotmatrix(V,F)
  % L = cotmatrix(V,T)
  %
  % For size(F,2)==4, This is distinctly NOT following definition that appears
  % in the appendix of: ``Interactive Topology-aware Surface Reconstruction,''
  % by Sharf, A. et al
  % http://www.cs.bgu.ac.il/~asharf/Projects/InSuRe/Insure_siggraph_final.pdf
  %
  % Instead it is a purely geometric construction. Find more details in Section
  % 1.1 of "Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
  % shapes" [Jacobson 2013]
  %
  % ND derivation given in "A MONOTONE FINITE ELEMENT SCHEME FOR
  % CONVECTION-DIFFUSION EQUATIONS" [Xu & ZIKATANOV 1999]
  %
  % 3D derivation given in "Aspects of unstructured grids and finite-volume
  % solvers for the Euler and Navier-Stokes equations" [Barth 1992]
  %
  %
  % Inputs:
  %   V  #V x dim matrix of vertex coordinates
  %   F  #F x simplex-size matrix of indices of triangle or tetrahedron corners
  % Outputs:
  %   L  sparse #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch), Denis Zorin
  %
  % See also: cotangent
  %

  ss = size(F,2);
  switch ss
  case 3
    %% Could just replace everything with:
    %C = cotangent(V,F);
    %L = sparse(F(:,[2 3 1]), F(:,[3 1 2]), C,size(V,1),size(V,1));
    %L = L+L';
    %L = L-diag(sum(L,2));
    
    % should change code below, so we don't need this transpose
    if(size(F,1) == 3)
      warning('F seems to be 3 by #F, it should be #F by 3');
    end
    F = F';

    % renaming indices of vertices of triangles for convenience
    i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
    % #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);
    % computing *unsigned* areas 
    if size(V,2) == 2
        % 2d vertex data
        dblA = abs(v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1));
    elseif size(V,2) == 3
        %n  = cross(v1,v2,2);  dblA  = multinorm(n,2);

        % area of parallelogram is twice area of triangle
        % area of parallelogram is || v1 x v2 || 
        n  = cross(v1,v2,2); 
        % THIS DOES MATRIX NORM!!! don't use it!!
        % dblA  = norm(n,2);

        % This does correct l2 norm of rows
        dblA = (sqrt(sum((n').^2)))';
    else 
        error('unsupported vertex dimension %d', size(V,2))
    end
    % cotangents and diagonal entries for element matrices
    cot12 = -dot(v1,v2,2)./dblA/2; cot23 = -dot(v2,v3,2)./dblA/2; cot31 = -dot(v3,v1,2)./dblA/2;
    % diag entries computed from the condition that rows of the matrix sum up to 1
    % (follows from  the element matrix formula E_{ij} = (v_i dot v_j)/4/A )
    diag1 = -cot12-cot31; diag2 = -cot12-cot23; diag3 = -cot31-cot23;
    % indices of nonzero elements in the matrix for sparse() constructor
    i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
    j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
    % values corresponding to pairs form (i,j)
    v = [cot12 cot12 cot23 cot23 cot31 cot31 diag1 diag2 diag3];
    % for repeated indices (i,j) sparse automatically sums up elements, as we
    % want
    L = sparse(i,j,v,size(V,1),size(V,1));
  case 4
    if(size(F,1) == 4 && size(F,2) ~=4)
      warning('F seems to be 4 by #F, it should be #F by 4');
    end
    % number of mesh vertices
    n = size(V,1);
    % cotangents of dihedral angles
    C = cotangent(V,F);
    %% TODO: fix cotangent to have better accuracy so this isn't necessary
    %% Zero-out almost zeros to help sparsity
    %C(abs(C)<10*eps) = 0;
    % add to entries
    L = sparse(F(:,[2 3 1 4 4 4]),F(:,[3 1 2 1 2 3]),C,n,n);
    % add in other direction
    L = L + L';
    % diagonal is minus sum of offdiagonal entries
    L = L - diag(sum(L,2));
    %% divide by factor so that regular grid laplacian matches finite-difference
    %% laplacian in interior
    %L = L./(4+2/3*sqrt(3));
    %% multiply by factor so that matches legacy laplacian in sign and
    %% "off-by-factor-of-two-ness"
    %L = L*0.5;
    % flip sign to match cotmatix.m
    if(all(diag(L)>0))
      warning('Flipping sign of cotmatrix3, so that diag is negative');
      L = -L;
    end
  end
end
