function [VV,FF] = remove_degenerate_faces(V,F,varargin)
  % REMOVE_DEGENERATE_FACES Remove degenerate faces from a mesh (V,F)
  % this can make combinatorially manifold meshes into non-manifold ones, but
  % should not change the "carrier" solid of solid meshes.
  % 
  % [VV,FF] = remove_degenerate_faces(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into V
  %   Optional:
  %      'MaxIter'  followed by the maximum number of iterations {100}
  %      'Epsilon'  followed by minimum area to be considered degenerate {0}
  % Outputs:
  %   VV  #VV by 3 list of vertex positions
  %   FF  #FF by 3 list of face indices into VV
  %

  % default values
  max_iter = 100;
  epsilon = 0;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','Epsilon'}, ...
    {'max_iter','epsilon'});
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

  %iter=1;
  %while true
  %  %% Snap exactly duplicate vertices
  %  %[V,~,J] = remove_duplicate_vertices(V,0);
  %  %F = J(F);
  %  %% Remove combinatorially degenerate faces
  %  %F = F(F(:,1)~=F(:,2) & F(:,2)~=F(:,3) & F(:,3)~=F(:,1),:);
  %  %% Remove geometrically degenerate faces
  %  %D = doublearea(V,F)<=epsilon;
  %  %if ~any(D)
  %  %  break;
  %  %end
  %  %% Find degenerate faces: since we've removed duplicate vertices these
  %  %% should only be obtuse triangles
  %  %FD = F(D,:);
  %  %% Flip longest edge (opposite obtuse angle)
  %  %allE = [FD(:,[2 3]);FD(:,[3 1]);FD(:,[1 2])];
  %  %l = edge_lengths(V,FD);
  %  %[~,longest] = max(l,[],2);
  %  %to_flip = allE((1:size(FD,1))'+size(FD,1)*(longest-1),:);
  %  %to_flip = unique(sort(allE((1:size(FD,1))'+size(FD,1)*(longest-1),:),2),'rows');
  %  %F = flip_edges(F,to_flip);
  %end



  % Force boundary vertices to come first so that snapping small edges to
  % minimum index will encourage snapping onto the boundary
  [V,FF,~,J] = faces_first(V,F,boundary_faces(F));
  EJ = (1:size(V,1))';

  % Combinatorially degenerate faces
  combinatorially_degenerate = ...
    @(G) any(G(:,2)==G(:,3) | G(:,3)==G(:,1) | G(:,1)==G(:,2),2);
  iter=1;
  get_edge = @(G,I,S) ...
      G([ ...
        sub2ind(size(G),find(I),mod(S+1-1,3)+1) ...
        sub2ind(size(G),find(I),mod(S+2-1,3)+1)]);
  validate_edges =  @(F,E) assert( ...
    all(ismember(sort(E,2),sort([F(:,[2 3]);F(:,[3 1]);F(:,[1 2])],2),'rows')));
  while true
    [IA,~,l,A] = is_acute(V,FF);
    D = A/2 < epsilon;
    % Find non-delaunay, longest edges of small, non-acute triangles
    L = ~is_delaunay(V,FF) & bsxfun(@eq,l,max(l,[],2)) & repmat(D&~IA,1,3);
    [LI,LJ] = ind2sub(size(FF),find(L));
    EL = FF([ ...
      sub2ind(size(FF),LI,mod(LJ+1-1,3)+1) ...
      sub2ind(size(FF),LI,mod(LJ+2-1,3)+1) ]);
    % Include all "non flippable" obtuse cases in collapsable
    IC = D & ~any(L,2);
    % Find short edges of acute triangles
    [~,S] = min(l(D&IC,:),[],2);
    ES = get_edge(FF,D&IC,S);

    if size(ES,1) == 0 && size(EL,1) == 0
      break;
    end

    validate_edges(FF,ES);
    validate_edges(FF,EL);

    %clf;
    %hold on;
    %tsurf(FF,V,'VertexIndices',1,'CData',1*D,'EdgeColor',[0.5 0.5 0.5]);
    %plot_edges(V,[ES;EL],'k','LineWidth',5);
    %plot_edges(V,ES,'r','LineWidth',3);
    %plot_edges(V,EL,'b','LineWidth',3);
    %hold off;
    %axis equal;
    %pause

    [ESL,EI] = conservative_edge_matching([ES;EL],'Method','recursive');
    %fprintf('%d+%d = %d %d\n',[size(ES,1) size(EL,1) size(ES,1)+size(EL,1) size(ESL,1)]);

    EL = ESL(EI>size(ES,1),:);
    ES = ESL(EI<=size(ES,1),:);

    validate_edges(FF,ES);
    validate_edges(FF,EL);


    %clf;
    %hold on;
    %tsurf(FF,V,'VertexIndices',1,'CData',1*D,'EdgeColor',[0.5 0.5 0.5]);
    %plot_edges(V,ESL,'k','LineWidth',5);
    %plot_edges(V,ES,'r','LineWidth',3);
    %plot_edges(V,EL,'b','LineWidth',3);
    %hold off;
    %axis equal;
    %pause

    % Flipping an edge can never create a degenerate face so do those first
    validate_edges(FF,EL);
    FF = flip_edges(FF,EL);
    validate_edges(FF,ES);
    ES = sort(ES,2);
    EJ(ES(:,2)) = ES(:,1);
    FF = EJ(FF);
    FF = FF(~combinatorially_degenerate(FF),:);

    if iter == max_iter
      warning('Maximum iterations reached');
      break;
    end
    iter = iter + 1;

  end

  [VV,UJ,I] = remove_unreferenced(V,FF);
  FF = UJ(FF);
  J = UJ(EJ(J));

  %clf;
  %hold on;
  %tsurf(FF,VV,'VertexIndices',1,'CData',1*D,'EdgeColor',[0.5 0.5 0.5]);
  %hold off;
  %axis equal;
  %pause
end
