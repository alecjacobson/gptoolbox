function M = massmatrix3(V,T, type)
  % MASSMATRIX3 mass matrix for the mesh given by V and F
  %
  % M = massmatrix3(V,T, type)
  %
  %
  % Inputs:
  %   V #V x 3 matrix of vertex coordinates
  %   T #T x 4  matrix of indices of tetrahedral corners
  %   type  string containing type of mass matrix to compute
  %     'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  %       of volumes of surrounding tets
  %     Not yet supported:
  %       'full': full mass matrix for p.w. linear fem
  %       'voronoi'
  % Output:
  %   M  #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: massmatrix
  %

  % vertices must be defined in 3D
  assert(size(V,2)==3);

  % should change code below, so we don't need this transpose
  if(size(T,1) == 4)
    warning('T seems to be 4 by #T, it should be #T by 4');
  end

  if strcmp(type,'full')
    error('full not supported yet...')
  elseif strcmp(type,'barycentric')
    %a = V(T(:,1),:);
    %b = V(T(:,2),:);
    %c = V(T(:,3),:);
    %d = V(T(:,4),:);
    %% http://en.wikipedia.org/wiki/Tetrahedron#Volume
    %% volume for each tetrahedron
    %v = repmat(abs(dot((a-d),cross2(b-d,c-d),2))./6./4,1,4);
    v = repmat(volume(V,T),1,4);

    % only diagonal elements
    i = T;
    M = sparse(T(:),T(:),v,size(V,1),size(V,1));
  elseif strcmp(type,'voronoi')
    pa = V(T(:,1),:);
    pb = V(T(:,2),:);
    pc = V(T(:,3),:);
    pd = V(T(:,4),:);
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
    M = sparse(T,T,pv,size(V,1),size(V,1));
  else 
    error('bad mass matrix type')
  end

  function r = cross2(a,b)
    % Optimizes r = cross(a,b,2), that is it computes cross products per row
    % Faster than cross if I know that I'm calling it correctly
    r =[a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
        a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
        a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  end

end
