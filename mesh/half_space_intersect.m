function [VV,FF,birth,UT,E] = half_space_intersect(V,F,p,n,varargin)
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
    [CV,CF] = cube(2,2,2);
    R = vrrotvec2mat(vrrotvec(n,[0 0 1]));
    % project vertex mean onto plane
    pp = mean(V) - ((mean(V)-p)*n')*n;
    % rotate box around projected mean
    CV = bsxfun(@plus,bsxfun(@minus,CV,[0.5 0.5 0])*diag([2 2 1])*bbd*R,pp);

    [VV,FF,J] = mesh_boolean(V,F,CV,CF,'intersect');
    if ~cap
      FF = FF(J<=size(F,1),:);
      [VV,IM] = remove_unreferenced(VV,FF);
      FF = IM(FF);
    elseif remesh_cap
      VA = VV;
      FA = FF(J<=size(F,1),:);
      JA = J(J<=size(F,1));
      VB = VV;
      FB = FF(J>size(F,1),:);
      [VB,FB] = remesh_planar_patches( ...
        VB,FB,'MinSize',4,'Force',true,'TriangleFlags','-q32Y');
      VV = [VA;VB];
      FF = [FA;size(VA,1)+FB];
      J = [JA;size(FA,1)+ones(size(FB,1),1)];
      [VV,~,IM] = remove_duplicate_vertices(VV,0);
      FF = IM(FF);
    end
    birth = J;
    birth(birth>size(F,1)) = 0;
  case 'fast'
    plane = [n(:)' -dot(p,n)];

    % Homogeneous coordinates
    IV = sum(bsxfun(@times,[V(:,1:3) ones(size(V,1),1)],plane),2);
    IF = IV(F);
    IF = reshape(IF,size(F));
    I13 = sum(IF<0,2) == 1;
    U = V;
    [U,E13,F13below,F13above,BC13] = one_below(U,F(I13,:),IF(I13,:));

    I31 = sum(IF>0,2) == 1;
    [U,E31,F31below,F31above,BC31] = one_below(U,F(I31,:),-IF(I31,:));

    VV = U;
    above = all(IF>=0,2);
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
      z = [0 0 1];
      c = dot(un,z);
      if abs(1-c)>eps
        v = cross(un,z);
        s = norm(v);
        cp = @(v) [0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
        R = eye(3) + cp(v) + cp(v)*cp(v)*(1-c)/s^2;
        T = [R [0;0;R(3,:)*(un'/norm(plane(1:3))*plane(4))]];
        UT = U(:,1:3)*T(1:2,1:3)';
      else
        UT = U(:,1:2);
      end

      if any(sum(adjacency_matrix(E),2) < 2)
        warning('Open boundary... cap will be wrong...');
      end
      [~,SF] = triangle(UT,E,[],'Flags','-Yc');
      RIM = full(sparse(IM,1,1:numel(IM)));
      % Index back onto U
      SF = RIM(SF);
      BC = barycenter(VV,SF);
      w = winding_number(V,F,BC);
      SF = SF(w>0.5,:);
      FF = [FF;fliplr(SF)];
      birth = [birth;zeros(size(SF,1),1)];
    end
    [VV,IM,J] = remove_unreferenced(VV,FF);
    FF = IM(FF);
  end

end
