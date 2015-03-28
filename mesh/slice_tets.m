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
  %   Optional:
  %     'Manifold' followed by whether to stitch together triangles into a
  %       manifold mesh {true}: results in more compact U but slightly slower.
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
  %     [U,G,BC] = slice_tets(V,T,[0 0 1 -z]);
  %     set(s,'Vertices',U,'Faces',G,'CData',BC*H);
  %     p.ZData = z*[1;1;1;1;1];
  %     drawnow;
  %   end
  %

  function [U,G,BC] = one_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    lambda = sIT(:,2:4)./bsxfun(@minus,sIT(:,2:4),sIT(:,1));
    if construct_BC
      BC = sparse( ...
        repmat((1:size(sT,1)*3)',1,2), ...
        [repmat(sT(:,1),3,1) reshape(sT(:,2:4),size(sT,1)*3,1)], ...
        [lambda(:) 1-lambda(:)], ...
        size(sT,1)*3,size(V,1));
      U = BC * V;
    else
      U = bsxfun(@times,V(repmat(sT(:,1),3,1),:),lambda(:)) + ...
        bsxfun(@times,1-lambda(:),V(sT(:,2:4),:));
      BC = [];
    end
    G = bsxfun(@plus,1:size(sT,1),[0;1;2]*size(sT,1))';
  end

  function [U,G,BC] = two_below(V,T,IT)
    [sIT,sJ] = sort(IT,2);
    sT = T(sub2ind(size(T),repmat(1:size(T,1),size(T,2),1)',sJ));
    lambda = sIT(:,3:4)./bsxfun(@minus,sIT(:,3:4),sIT(:,1));
    gamma  = sIT(:,3:4)./bsxfun(@minus,sIT(:,3:4),sIT(:,2));
    if construct_BC
      BC = sparse( ...
        repmat((1:size(sT,1)*4)',1,2), ...
        [repmat(sT(:,1),2,1) reshape(sT(:,3:4),size(sT,1)*2,1); ...
         repmat(sT(:,2),2,1) reshape(sT(:,3:4),size(sT,1)*2,1)], ...
        [lambda(:) 1-lambda(:);gamma(:) 1-gamma(:)], ...
        size(sT,1)*4,size(V,1));
      U = BC * V;
    else 
      U = [ ...
        bsxfun(@times,V(repmat(sT(:,1),2,1),:),lambda(:)) + ...
        bsxfun(@times,1-lambda(:),V(sT(:,3:4),:)); ...
        bsxfun(@times,V(repmat(sT(:,2),2,1),:),gamma(:)) + ...
        bsxfun(@times,1-gamma(:),V(sT(:,3:4),:))];
      BC = [];
    end
    G = [ ...
      bsxfun(@plus,1:size(sT,1),[0;1;3]*size(sT,1))'; ...
      bsxfun(@plus,1:size(sT,1),[0;3;2]*size(sT,1))'];
  end

  % default values
  manifold = true;
  construct_BC = nargout >= 4;

  % Map of parameter names to variable names
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
  IV = sum(bsxfun(@times,[V ones(size(V,1),1)],plane),2);
  IT = IV(T);
  IT = reshape(IT,size(T));

  I13 = sum(IT<0,2) == 1;
  [U13,G13,BC13] = one_below(V,T(I13,:),IT(I13,:));
  I31 = sum(IT>0,2) == 1;
  [U31,G31,BC31] = one_below(V,T(I31,:),-IT(I31,:));
  I22 = sum(IT<0,2) == 2;
  [U22,G22,BC22] = two_below(V,T(I22,:),IT(I22,:));

  U = [U13;U31;U22];
  if construct_BC
    BC = [BC13;BC31;BC22];
  end
  G = [G13;size(U13,1)+[G31;size(U31,1)+[G22;]]];
  J = [find(I13);find(I31);repmat(find(I22),2,1)];
  N = normals(U,G);
  flip = sum(bsxfun(@times,N,plane(1:3)),2)<0;
  G(flip,:) = fliplr(G(flip,:));

  if manifold
    % should be able to do this combinatorially
    bbd = normrow(max(V)-min(V));
    [U,I,IM] = remove_duplicate_vertices(U,1e-14*bbd);
    if construct_BC
      BC = BC(I,:);
    end
    G = IM(G);
  end

end
