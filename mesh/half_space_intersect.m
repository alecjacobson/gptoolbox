function [VV,FF,birth] = half_space_intersect(V,F,p,n,varargin)
  % HALF_SPACE_INTERSECT Intersect a closed mesh (V,F) with a half space
  % defined by a point on a plane the plane's normal.
  %
  % [VV,FF] = half_space_intersect(V,F,p,n)
  % [VV,FF] = half_space_intersect(V,F,p,n,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  %   p  1 by 3 point on plane position
  %   n  1 by 3 plane normal vector
  %   Optional:
  %     'Cap' followed by whether to include a cap along plane {true}
  %     'Manifold' followed by whether to stitch together triangles into a
  %       manifold mesh {true}: results in more compact U but slightly slower.
  %     'Method' Followed by either:
  %        {'fast'}  quickly mesh all triangles crossing plane (computed via
  %          inexact floating point operations).
  %        'exact'  literally construct a plane and intersect the input mesh
  %          against it.
  %     'RemeshCap' followed by wether to remesh cap for better triangulation
  %       (only makes sense for 'Method','exact') {false}
  % Outputs:
  %   VV  #VV by 3 list of new mesh vertex positions
  %   FF  #FF by 3 list of new mesh triangle indices
  %   birth   #FF list of indices into F of birth parents, 0 means on cap
  %

  function [U,E,Fbelow,Fabove,BC] = one_below(V,F,IF)
    [~,J] = min(IF,[],2);
    I = sub2ind(size(F),repmat(1:size(F,1),size(F,2),1)',mod([J J+1 J+2]-1,3)+1);
    lambda = IF(I(:,2:3))./bsxfun(@minus,IF(I(:,2:3)),IF(I(:,1)));
    BC = sparse( ...
      repmat((1:size(F,1)*2)',1,2), ...
      [repmat(F(I(:,1)),2,1) reshape(F(I(:,2:3)),size(F,1)*2,1)], ...
      [lambda(:) 1-lambda(:)], ...
      size(F,1)*2,size(V,1));
    E = size(V,1)+bsxfun(@plus,1:size(F,1),[0;1]*size(F,1))';
    U = [V;BC * V];
    Fbelow = [F(I(:,1)) E];
    Fabove = [F(I(:,2)) fliplr(E);F(I(:,2:3)) E(:,2)];
  end

  method = 'fast';
  remesh_cap = false;
  construct_BC = false;
  cap = true;
  manifold = true;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Manifold','Method','RemeshCap','Cap'}, ...
    {'manifold','method','remesh_cap','cap'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  switch method
  case 'exact'
    bbd = sqrt(sum((max(V)-min(V)).^2,2));
    % row vectors
    n = reshape(n,1,[]);
    p = reshape(p,1,[]);
    N = null(n)';
    CV = bsxfun(@plus,p,bbd*[1 1;1 -1;-1 -1;-1 1]*N);
    CF = [1 2 3;1 3 4];
  
    VCV = [V;CV];
    FCF = [F;size(V,1)+CF];
    [VV,FF,~,J,IM] = selfintersect(VCV,FCF);
  
    if cap
      CF = FF(J>size(F,1),:);
      BC = barycenter(VV,CF);
      w = winding_number(V,F,BC);
      CF = CF(w>0.5,:);
    else
      CF = [];
    end
  
    F = FF(J<=size(F,1),:);
    birth = J(J<=size(F,1));
    BC = barycenter(VV,F);
    above = sum(bsxfun(@times,bsxfun(@minus,BC,p),n),2)>=0;
    F = F(above,:);
    birth = birth(above);
    FF = [F;CF];
    birth = [birth;zeros(size(CF,1),1)];
    [VV,~,J] = remove_duplicate_vertices(VV,1e-10);
    FF = J(FF);
    [VV,J] = remove_unreferenced(VV,FF);
    FF = J(FF);
    if remesh_cap && cap
      [VV,FF] = remesh_planar_patches(VV,FF,'MinSize',4,'Force',true, ...
        'Except',1:size(F,1));
      % better have only remeshed on cap.
      birth = [birth(1:size(F,1));zeros(size(FF,1)-size(F,1),1)];
    end
  case 'fast'
    plane = [n(:)' -dot(p,n)];

    % Homogeneous coordinates
    IV = sum(bsxfun(@times,[V ones(size(V,1),1)],plane),2);
    IF = IV(F);
    IF = reshape(IF,size(F));
    I13 = sum(IF<0,2) == 1;
    U = V;
    [U,E13,F13below,F13above,BC13] = one_below(U,F(I13,:),IF(I13,:));

    I31 = sum(IF>0,2) == 1;
    [U,E31,F31below,F31above,BC31] = one_below(U,F(I31,:),-IF(I31,:));

    VV = U;
    above = all(IF>0,2);
    FF = [F(above,:);F13above;F31below];
    birth = [find(above);repmat(find(I13),2,1);find(I31)];
    E = [E13;E31];

    if manifold
      % should be able to do this combinatorially
      bbd = normrow(max(V)-min(V));
      [VV,I,IM] = remove_duplicate_vertices(VV,1e-14*bbd);
      if construct_BC
        BC = BC(I,:);
      end
      FF = IM(FF);
      E = IM(E);
    end

    if cap
      [U,IM] = remove_unreferenced(VV,E);
      E = IM(E);
      % mesh of cap
      % compute transformation taking plane to (XY)-plane
      un  = plane(1:3)/norm(plane(1:3));
      v = cross(un,[0 0 1]);
      c = dot(un,[0 0 1]);
      s = norm(v);
      cp = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
      R = eye(3) + cp(v) + cp(v)*cp(v)*(1-c)/s^2;
      T = [R [0;0;R(3,:)*(un'/norm(plane(1:3))*plane(4))]];
      [~,SF] = triangle(U*T(1:2,1:3)',E,[],'Flags','-Y');
      RIM = full(sparse(IM,1,1:numel(IM)));
      % Index back onto U
      SF = RIM(SF);
      BC = barycenter(VV,SF);
      w = winding_number(V,F,BC);
      SF = SF(w>0.5,:);
      FF = [FF;fliplr(SF)];
      birth = [birth;zeros(size(SF,1),1)];
    end
    [VV,IM] = remove_unreferenced(VV,FF);
    FF = IM(FF);
  end

end
