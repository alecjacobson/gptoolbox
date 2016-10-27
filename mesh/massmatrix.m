function M = massmatrix(V,F, type)
  % MASSMATRIX mass matrix for the mesh given by V and F
  %
  % M = massmatrix(V,F, type)
  % M = massmatrix(V,T, type)
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x simplex-size  matrix of indices of triangle corners
  %  type  string containing type of mass matrix to compute
  %   'full': full mass matrix for p.w. linear fem
  %   'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  %   'voronoi': true voronoi area, except in cases where triangle is obtuse
  %     then uses 1/2, 1/4, 1/4 {simplex size 3 only}
  % Output:
  %  M  #V by #V sparse mass matrix
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  function r = cross2(a,b)
    % Optimizes r = cross(a,b,2), that is it computes cross products per row
    % Faster than cross if I know that I'm calling it correctly
    r =[a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
        a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
        a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  end

  % simplex size
  ss = size(F,2);
  if nargin<3
    switch ss
    case 3 
      type = 'voronoi';
    case 4
      type = 'barycentric';
    end
  end

  switch ss
  case 3
    % should change code below, so we don't need this transpose
    if(size(F,1) == 3)
      warning('F seems to be 3 by #F, it should be #F by 3');
    end
    F = F';
    
    % renaming indices of vertices of triangles for convenience
    i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
    % #F x 3 matrices of triangle edge vectors, named after opposite vertices
    v1 = V(i3,:) - V(i2,:);  v2 = V(i1,:) - V(i3,:); v3 = V(i2,:) - V(i1,:);
    % computing the areas
    if size(V,2) == 2
    % 2d vertex data
      dblA = v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1);
    elseif size(V,2) == 3
      %n  = cross(v1,v2,2);  dblA  = multinorm(n,2);
      n  = cross(v1,v2,2);  
      
      % dblA  = norm(n,2);
      % This does correct l2 norm of rows
      dblA = (sqrt(sum((n').^2)))';
    else 
      error('unsupported vertex dimension %d', size(V,2))
    end
    if strcmp(type,'full')
        % arrays for matrix assembly using 'sparse'
        % indices and values of the element mass matrix entries in the order 
        % (1,2), (2,1),(2,3), (3,2), (3,1), (1,3) (1,1), (2,2), (3,3);
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
        M = sparse(i,j,v,size(V,1), size(V,1));
    elseif strcmp(type,'barycentric')
        % only diagonal elements
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v,diag_v,diag_v];
        M = sparse(i,j,v,size(V,1), size(V,1));
    elseif strcmp(type,'voronoi')

      % just ported version of intrinsic code

      % edges numbered same as opposite vertices
      FT = F';
      l = [ ...
        sqrt(sum((V(FT(:,2),:)-V(FT(:,3),:)).^2,2)) ...
        sqrt(sum((V(FT(:,3),:)-V(FT(:,1),:)).^2,2)) ...
        sqrt(sum((V(FT(:,1),:)-V(FT(:,2),:)).^2,2)) ...
        ];
      M = massmatrix_intrinsic(l,F',size(V,1),'voronoi');
    else 
        error('bad mass matrix type')
    end
    
    % warn if any rows are all zero (probably unreferenced vertices)
    if(any(sum(M,2) == 0))
      warning('Some rows have all zeros... probably unreferenced vertices..');
    end
  case 4
    % vertices must be defined in 3D
    assert(size(V,2)==3);
  
    % should change code below, so we don't need this transpose
    if(size(F,1) == 4)
      warning('F seems to be 4 by #F, it should be #F by 4');
    end
  
    switch type
    case 'full'
      % e.g., from _Finite Elements for Analysis and Design_,  J. E. Akin, pp
      % 198. Or Zienkiewicz-4
      %
      % Earliest reference I could find is "Numerical Integration over
      % Simplexes and Cones", [Hammer et al. 1956]
      % 
      I = [F(:,[2 3 4 3 4 1 4 1 2 1 2 3 1 2 3 4])];
      J = [F(:,[1 1 1 2 2 2 3 3 3 4 4 4 1 2 3 4])];
      vol = abs(volume(V,F));
      VV = [repmat(vol/20,1,3*4) repmat(vol/10,1,4)];
      M = sparse(I,J,VV,size(V,1),size(V,1));
      % Uniform mid-face quadrature rules don't seem to be quadratically
      % precise
      % I = [F(:,[2 3 4 3 4 1 4 1 2 1 2 3 1 2 3 4])];
      % J = [F(:,[1 1 1 2 2 2 3 3 3 4 4 4 1 2 3 4])];
      % vol = abs(volume(V,F));
      % VV = [repmat(vol/18,1,3*4) repmat(vol/12,1,4)];
      % M = sparse(I,J,VV,size(V,1),size(V,1));
    case 'barycentric'
      %a = V(F(:,1),:);
      %b = V(F(:,2),:);
      %c = V(F(:,3),:);
      %d = V(F(:,4),:);
      %% http://en.wikipedia.org/wiki/Tetrahedron#Volume
      %% volume for each tetrahedron
      %v = repmat(abs(dot((a-d),cross2(b-d,c-d),2))./6./4,1,4);
      v = repmat(abs(volume(V,F)),1,4);
  
      % only diagonal elements
      M = sparse(F(:),F(:),v/4,size(V,1),size(V,1));
    case 'voronoi'
      pa = V(F(:,1),:);
      pb = V(F(:,2),:);
      pc = V(F(:,3),:);
      pd = V(F(:,4),:);
      % as if pd is the origin
      a = pa-pd;
      b = pb-pd;
      c = pc-pd;
      % circumcenter:
      % http://en.wikipedia.org/wiki/Tetrahedron#More_vector_formulas_in_a_general_tetrahedron
      cc = pd + bsxfun(@rdivide, ...
        bsxfun(@times,sum(a.*a,2),cross2(b,c)) + ...
        bsxfun(@times,sum(b.*b,2),cross2(c,a)) + ...
        bsxfun(@times,sum(c.*c,2),cross2(a,b)), ...
        2*sum(a.*cross2(b,c),2));
      % get correct sign
      sa = sign(sum((pa-pb).*cross2(pc-pb,pd-pb),2)/6);
      sb = sign(sum((pb-pc).*cross2(pd-pc,pa-pc),2)/6);
      sc = sign(sum((pc-pd).*cross2(pa-pd,pb-pd),2)/6);
      sd = sign(sum((pd-pa).*cross2(pb-pa,pc-pa),2)/6);
      % get area of sub-tet and correct sign
      la = sa.*sum((cc-pb).*cross2(pc-pb,pd-pb),2)/6;
      lb = sb.*sum((cc-pc).*cross2(pd-pc,pa-pc),2)/6;
      lc = sc.*sum((cc-pd).*cross2(pa-pd,pb-pd),2)/6;
      ld = sd.*sum((cc-pa).*cross2(pb-pa,pc-pa),2)/6;
      v = abs(sum(a.*cross2(b,c),2)/6);
      %max(abs((la+lb+lc+ld - v)))
      assert(all((la+lb+lc+ld - v)<10*eps));
      % partial volumes attached to each corner
      pv = [(lb+lc+ld)/3 (lc+ld+la)/3 (ld+la+lb)/3 (la+lb+lc)/3];
      M = sparse(F,F,pv,size(V,1),size(V,1));
    otherwise
      error('bad mass matrix type')
    end
  
  end
end
