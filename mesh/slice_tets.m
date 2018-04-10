function [U,G,J,BC] = slice_tets(V,T,plane,varargin)
  % SLICE_TETS Slice through a tet mesh (V,T) along a given plane (via its
  % implicit equation).
  %
  % [U,G] = slice_tets(V,T,plane)
  % [U,G,J,BC] = slice_tets(V,T,plane,'ParameterName',parameter_value, ...)
  %
  % Inputs:
  %   V  #V by 3 list of tet mesh vertices
  %   T  #T by 4 list of tet indices into V 
  %   plane  list of 4 coefficients in the plane equation: [x y z 1]'*plane = 0
  %     or
  %   S  #V by 1 list of values per vertex
  % Outputs:
  %   U  #U by 3 list of triangle mesh vertices along slice
  %   G  #G by 3 list of triangles indices into U
  %   J  #G list of indices into T revealing which tet this face came from
  %   BC  #U by #V list of barycentric coordinates (or more generally: linear
  %     interpolation coordinates) so that U = BC*V
  % 
  % Example:
  %   % Tet mesh in (V,T)
  %   F = boundary_faces(T);
  %   % Solve poisson equation for interesting function inside
  %   L = cotmatrix(V,T);
  %   M = massmatrix(V,T);
  %   b = unique(F);
  %   int = setdiff(1:size(V,1),b);
  %   H = zeros(size(V,1),1);
  %   H(int) = (-L(int,int))\(M(int,int)*ones(numel(int),1));
  %   clf;
  %   t = tsurf(F,V,'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.2,'EdgeAlpha',0.2);
  %   hold on;
  %     s = tsurf([1 1 1],V,'EdgeColor','none',fphong);
  %     BB = bounding_box(V(:,1:2));
  %     BB = BB([1 2 4 3 1],:);
  %     p = plot3(BB(:,1),BB(:,2),min(V(:,3))*[1;1;1;1;1],'-','LineWidth',3);
  %   hold off;
  %   caxis([min(H) max(H)]);
  %   axis equal;
  %   for z = linspace(min(V(:,3)),max(V(:,3)))
  %     [U,G,J,BC] = slice_tets(V,T,[0 0 1 -z]);
  %     set(s,'Vertices',U,'Faces',G,'CData',BC*H);
  %     p.ZData = z*[1;1;1;1;1];
  %     drawnow;
  %   end
  %
  % Example:
  %   % slice a triangle mesh to produce an oriented curve
  %   plane = [0 0 1 -40];
  %   [U,E,J] = slice_tets(V,F(:,[1 2 3 3]),plane);
  %   [U,~,I] = remove_duplicate_vertices(U,eps);
  %   E(E(:,2) == E(:,1),2) = E(E(:,2) == E(:,1),3);
  %   E = E(:,1:2);
  %   N = normals(V,F(J,:));
  %   N = N-sum(N.*plane(1:3),2).*plane(1:3);
  %   M = U(E(:,2),:)-U(E(:,1),:);
  %   Q = [1 plane(1:3)].*[cos(pi/4) sin(pi/4)*[1 1 1]];
  %   W = quatmultiply(quatmultiply(Q,[zeros(size(M,1),1) M]),Q.*[1 -1 -1 -1]);
  %   W = W(:,2:4);
  %   R = sign(sum(W.*N,2))>0;
  %   E(R,:) = fliplr(E(R,:));
  %

  flipped_order = flipped_tet_orders();

  function [U,G] = one_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    U = [repmat(sT(:,1),3,1) reshape(sT(:,2:4),size(sT,1)*3,1)];
    G = bsxfun(@plus,1:size(sT,1),[0;1;2]*size(sT,1))';
    flip = ismember(sJ,flipped_order,'rows');
    G(flip,:) = fliplr(G(flip,:));
  end

  function [U,G] = two_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    U =  ...
        [repmat(sT(:,1),2,1) reshape(sT(:,3:4),size(sT,1)*2,1); ...
         repmat(sT(:,2),2,1) reshape(sT(:,3:4),size(sT,1)*2,1)];
    G = [ ...
      bsxfun(@plus,1:size(sT,1),[0;1;3]*size(sT,1))'; ...
      bsxfun(@plus,1:size(sT,1),[0;3;2]*size(sT,1))'];
    flip = ismember([sJ;sJ],flipped_order,'rows');
    G(flip,:) = fliplr(G(flip,:));
  end

  % default values
  manifold = true;
  construct_BC = nargout >= 4;

  % Map of parameter names to variable names
  % 'Manifold' is kept for legacy reason, but is no longer needed (always
  % "manifold")
  params_to_variables = containers.Map( ...
    {'Manifold'}, ...
    {'manifold'});
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

  % Homogeneous coordinates
  if numel(plane) == 4
    IV = sum(bsxfun(@times,[V ones(size(V,1),1)],plane),2);
  else
    IV = plane;
  end
  IT = IV(T);
  IT = reshape(IT,size(T));

  I13 = sum(IT<0,2) == 1;
  [U13,G13] = one_below(V,T(I13,:),IT(I13,:));
  I31 = sum(IT>0,2) == 1;
  [U31,G31] = one_below(V,T(I31,:),-IT(I31,:));
  I22 = sum(IT<0,2) == 2;
  [U22,G22] = two_below(V,T(I22,:),IT(I22,:));

  U = [U13;U31;U22];

  % sort edges in U
  sU = sort(U,2);
  % find unique edges in U
  [E,uI,uJ] = unique(sU,'rows');
  % Values at each corner of each unique edge
  IE = IV(E);
  % Sort endpoints of each edge by value
  [sIE,sJ] = sort(IE,2);
  sE = E(sub2ind(size(E),repmat(1:size(E,1),size(E,2),1)',sJ));
  % Lambda from smallest to largest endpoint
  lambda = sIE(:,2)./(sIE(:,2)-sIE(:,1));
  % Vertex position on each unique edge
  U = V(sE(:,1),:).*lambda+ V(sE(:,2),:).*(1-lambda);
  G = [G13;size(U13,1)+[fliplr(G31);size(U31,1)+[G22;]]];
  G = uJ(G);

  if construct_BC
    BC = sparse( ...
      repmat(1:size(sE,1),2,1)', ...
      sE, ...
      [lambda 1-lambda], ...
      size(sE,1),size(V,1));
  end
  J = [find(I13);find(I31);repmat(find(I22),2,1)];

end
