function [V,F,b,bc] = triangulate_curves(P,varargin)
  % TRIANGULATE_CURVES Construct a (properly slitted) mesh for a given set of curves and
  % determine corresponding boundary conditions for solving diffusion curves.
  % 
  % [V,F,b,bc] = triangulate_curves(P,varargin)
  %
  % Inputs:
  %    #curve by 1 list of lists of curve points
  %  Optional:
  %    'BoundingBox' followed by whether to include bounding box (otherwise
  %    there better be an outer loop.
  %
  % Example:
  %     clf;
  %     % get curve (user draws)
  %     Praw = {};
  %     while true
  %       [Prawc,p] = get_pencil_curve();
  %       if size(Prawc,1) == 1 || sum((max(Prawc)-min(Prawc)).^2,2)<1e-5
  %         fprintf('continuing\n');
  %         break;
  %       end
  %       Praw{end+1} = Prawc;
  %       %% remove plot
  %       %delete(p);
  %     end
  %     % Simplify
  %     P = cell(numel(Praw),1);
  %     for c = 1:numel(Praw)
  %       P{c} = dpsimplify(Praw{c},0.001);
  %     end
  %     [V,F,b,bc] = triangulate_curves(P);
  %     W = harmonic(V,F,b,bc);
  %     colors = random_color(4*numel(P));
  %     RGB = W * colors;
  %     tsurf(F,[V W(:,4)],'FaceVertexCData',clamp(RGB),fphong,'EdgeColor','none');axis equal;view(2)
  %

  % default values
  use_bounding_box = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'BoundingBox'}, ...
    {'use_bounding_box'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  E = [];
  % list of all points
  PP = [];
  % list of curve sizes to determine from boundary markers which curve edge came
  % from
  curve_pos_ind = 0;
  for c = 1:numel(P)
    if isstruct(P{c})
      Pc = P{c}.P;
      Ec = P{c}.E;
    else
      Pc = P{c};
      Ec = [1:size(Pc,1)-1;2:size(Pc,1)]';
    end
    % inefficient appending
    curve_pos_ind = [curve_pos_ind curve_pos_ind(end)+size(Ec,1)];
    E = [E;size(PP,1)+Ec];
    PP = [PP;Pc];
  end
  l = normrow(PP(E(:,1),:)-PP(E(:,2),:));

  Facets = [];
  Facets.facets = mat2cell(E,ones(size(E,1),1),2);
  Facets.boundary_marker = (1:size(E,1))'+1;

  % need unique set of points
  assert(size(PP,1) == size(unique(PP,'rows'),1));
  max_area = min(l)^2;
  for pass = 1:2
    prefix = tempprefix();
    if use_bounding_box
      BB = bounding_box(PP);
      BB = bsxfun(@plus,1.1*bsxfun(@minus,BB,mean(BB)),mean(BB));
      bb_flag = 'c';
    else 
      BB = [];
      bb_flag = '';
    end
    writePOLY_triangle([prefix '.poly'],[PP;BB],Facets,[]);
    % Careful, triangle does not understand scientific notation
    flags = sprintf('-q33pa%0.17f%s',max_area,bb_flag);
    %flags = '-cp';
    command = [path_to_triangle ' ' flags ' ' prefix];
    [status, result] = system( command );
    if(status ~= 0)
       error(result);
    end
    [V,I] = readNODE([prefix '.1.node']);
    F = readELE([prefix '.1.ele']);
    max_area = mean(doublearea(V,F))*2;
  end
  [~,BE,BM] = readPOLY_triangle([prefix '.1.poly']);
  BE = BE(BM>1,:);
  BM = BM(BM>1,:)-1;
  % Re-orient edges to follow original curve
  BEdotE = sum((V(BE(:,1),:)-V(BE(:,2),:)).*(PP(E(BM,1),:)-PP(E(BM,2),:)),2);
  BE(BEdotE<=0,:) = fliplr(BE(BEdotE<=0,:));
  uBE = unique(sort(BE,2),'rows');

  % boundary conditions (to be thinned)
  % #V by #curves * 2 (left/right) * 2 (source-to-dest/dest-to-source)
  bc = sparse(size(V,1),(numel(curve_pos_ind)-1)*2*2);

  paths = [];
  for c = 1:numel(P)
    % Boundary edges of this curve
    BEp = BE(BM>curve_pos_ind(c) & BM<=curve_pos_ind(c+1),:);
    BMp = BM(BM>curve_pos_ind(c) & BM<=curve_pos_ind(c+1));
    % Assumes E and thus BM are in order 
    % source of curve
    sc = E(curve_pos_ind(c)+1,1);
    path = [sc];
    % loop over original curve
    for e = (curve_pos_ind(c)+1):curve_pos_ind(c+1)
      % dest of last original edge is source of next
      se = path(end);
      assert(se == E(e,1));
      de = E(e,2);
      BEe = BEp(BMp==e,:);
      A = adjacency_matrix(BEe);
      [~,pathe] = graphshortestpath(A,se,de);
      path = [path pathe(2:end)];
    end
    % Arc length parameterization of path
    path_edge_lens = sqrt(sum((V(path(2:end),:)-V(path(1:end-1),:)).^2,2));
    path_t = [0;cumsum(path_edge_lens)./sum(path_edge_lens)];

    bc(path,(c-1)*2*2+([1 3])) = [path_t path_t];
    bc(path,(c-1)*2*2+([2 4])) = 1-[path_t path_t];
    % inefficient append
    paths = [paths path];
  end

  % remove intersection points
  [upaths,~,IC] = unique(paths);
  counts = histc(paths,upaths);
  bc(upaths(counts>1),:) = 0;

  % disconnected mesh with each corner as distinct vertex
  W = V(F(:),:);
  m = size(F,1);
  G = bsxfun(@plus,[0 m 2*m],[1:m]');
  % remap boundary conditions to new disconnected mesh
  bc = bc(F(:),:);

  % boundary
  left = reshape( ...
    ismember([F(:,[2 3]);F(:,[3 1]);F(:,[1 2])],BE,'rows'),m,3);
  right = reshape( ...
    ismember([F(:,[2 3]);F(:,[3 1]);F(:,[1 2])],fliplr(BE),'rows'),m,3);
  [touch,touchloc] = ismember(F,BE(:));
  touch = reshape(touch,m,3);
  % definitely left or right
  b_left  = G(  left(:,[2 3 1])| left(:,[3 1 2]));
  b_right = G(right(:,[2 3 1])|right(:,[3 1 2]));
  b_left = b_left(:);
  b_right = b_right(:);


  %Z = repmat(1:size(G,1),1,3)';
  %tsurf(G,[W Z],'FaceColor','g');
  %hold on;plot(P(:,1),P(:,2),'-o','LineWidth',4); hold off;
  %hold on;plot_edges(V,BE,'--r','LineWidth',4); hold off;
  %hold on; plot3(W(b_right,1),W(b_right,2),Z(b_right,1),'or','LineWidth',4);hold off
  %hold on; plot3(W(b_left,1),W(b_left,2),Z(b_left,1),'ob','Linewidth',2);hold off
  %error


  % vertex to face incidence
  %V2F = sparse( ...
  %  F(:), ...
  %  repmat(1:size(F,1),1,3)', ...
  %  [ones(size(F,1),1); 2*ones(size(F,1),1); 3*ones(size(F,1),1)], ...
  %  size(V,1),size(F,1));
  % Finding on the transpose is faster
  V2FT = sparse( ...
    repmat(1:size(F,1),1,3)', ...
    F(:), ...
    [ones(size(F,1),1); 2*ones(size(F,1),1); 3*ones(size(F,1),1)], ...
    size(F,1), ...
    size(V,1));

  %J = 1:size(W,1);
  W2V = sparse((1:size(F,1)*3)',F(:),1,size(W,1),size(V,1));
  % http://stackoverflow.com/a/15719703/148668
  [~,first] = max((W2V*W2V')~=0,[],2);
  assert(numel(first)==size(W,1));
  assert(all(first));
  % default is to map to first occurrence
  J = first;

  VB = unique(BE)';
  % loop over original vertices on curves
  for v = VB
    % Faces incident on v and index of corner in face
    %[~,IFv,CFv] = find(V2F(v,:));
    [IFv,~,CFv] = find(V2FT(:,v));
    Fv = F(IFv,:);
    [~,uE2Fv,uE] = edge_adjacency_matrix(Fv);
    isBE = ismember(uE,uBE,'rows');
    A = uE2Fv(~isBE,:)' * uE2Fv(~isBE,:);
    [~,C] = conncomp(A);
    % loop over components
    for c = 1:max(C)
      min_fv = find(C==c,1);
      % first index in w for this component
      min_w = G(IFv(min_fv),CFv(min_fv));
      % map all intances in this component to first 
      J(G(sub2ind(size(G),IFv(C==c),CFv(C==c)))) = min_w;
    end
  end

  % remap all at once
  G = J(G);

  b_left = unique(J(b_left)');
  b_right = unique(J(b_right)');
  [W,I] = remove_unreferenced(W,G);
  % Remap boundary conditions
  bc(I,:) = bc;
  % Only keep referenced
  bc = bc(1:size(W,1),:);
  % Remap faces
  G = I(G);
  % remap left/right
  b_left = I(b_left);
  b_right = I(b_right);

  % Let shared vertices---internal curve endpoints---"float"
  % But really this should be left up to implementation
  b_shared = intersect(b_left,b_right); 
  b_left = setdiff(b_left,b_shared);
  b_right = setdiff(b_right,b_shared);
  % zap left, right and shared accordingly
  bc(b_right,[1:4:end 2:4:end]) = 0;
  bc(b_left,[3:4:end 4:4:end]) = 0;
  bc(b_shared,:) = 0;

  has_b = any(bc,2);
  b = find(has_b);
  bc = bc(b,:);

  % rename
  V = W;
  F = G;

end
