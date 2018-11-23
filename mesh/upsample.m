function [VV,FF,FO] = upsample(V,F,varargin)
  % UPSAMPLE Upsample a mesh by adding vertices on existing edges/faces After n
  % iterations of loop subivision, the resulting mesh will have
  %                                      4^n |F| faces
  %            2^n*|E| +  3*2^(n-1)*(2^n-1))*|F| edges
  %   |V| + (2^n-1)|E| + (1+2^(n-1)*(2^n-3))*|F| vertices
  %
  % [VV,FF] = upsample(V,F)
  % [VV,FF,FO] = upsample(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by simplex-size list of simplex indices
  %   Optional:
  %     'KeepDuplicates' followed by either true or {false}}
  %     'OnlySelected' followed by a list of simplex indices into F to
  %       subdivide.
  %     'Iterations' followed by number of recursive calls {1}
  % Outputs:
  %  VV  #VV by dim list new vertex positions, original V always comes first.
  %  FF  #FF by simplex-size new list of face indices into VV
  %  FO  #FF list of indices into F of original "parent" simplex
  %  
  % This is Loop subdivision without moving the points
  %
  % Example:
  %   [VVB,FF] = upsample([V speye(size(V,1))],F);
  %   VV = full(VVB(:,1:size(V,2)));
  %   B = VVB(:,size(V,2)+1:end);
  %   max(max(abs(VV - B*V)))
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  
%   % Add a new vertex at each face barycenter
%   % compute barycenters
%   C = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
%   % append barycenters to list of vertex positions
%   VV = [V;C];
%   % list of indices to barycenters
%   i = size(V,1) + (1:size(C,1))';
%   % New face indices, 3 new face for each original face
%   FF = [F(:,1) F(:,2) i; F(:,2) F(:,3) i; F(:,3) F(:,1) i];
%   

  %      o           o     
  %     / \         / \     
  %    x   x  ---> o---o    
  %   /     \     / \ / \  
  %  o---x---o   o---o---o  
  function [U14,F14,E14] = one_four(offset,V,F)
    % compute midpoints (actually repeats, one midpoint per edge per face)
    E14 = [F(:,2) F(:,3);F(:,3) F(:,1);F(:,1) F(:,2)];
    U14 = (V(E14(:,1),:)+V(E14(:,2),:))/2;
    % indices of midpoints
    nu = size(U14,1);
    i1 = offset+(1:(nu/3))';
    i2 = offset+((nu/3)) + (1:(nu/3))';
    i3 = offset+((nu/3)+(nu/3)) + (1:(nu/3))';
    % new face indices, 4 new faces for each original face. As if we simply
    % ignored the duplicates in m and had appended m to V
    F14 = [ F(:,1) i3 i2 ; F(:,2) i1 i3 ; F(:,3) i2 i1 ; i1 i2 i3];
  end

  %      o           o     
  %     / \         /|\     
  %    x   \  ---> o | \    
  %   /     \     / \|  \  
  %  o---x---o   o---o---o  
  function [U13,F13,E13] = one_three(offset,V,F,M)
    E = [F(:,2) F(:,3);F(:,3) F(:,1);F(:,1) F(:,2)];
    A = cumsum(M,2);
    [SJ,SI] = find(M'==0);
    % Vertex opposite non-subdivided edge
    flip = SJ==2;
    O1  = F(sub2ind(size(A),SI,SJ));
    % Next vertex
    O2  = F(sub2ind(size(A),SI,mod(SJ,3)+1));
    % Next next vertex
    O3  = F(sub2ind(size(A),SI,mod(SJ+1,3)+1));
    % Vertex opposite first subdivided edge
    first = M'==1 & A'==1;
    [SJ,SI] = find(first);
    I1 = sub2ind(size(A),SI,SJ);
    E13 = [ ...
      F(sub2ind(size(A),SI,mod(SJ,3)+1)) ...
      F(sub2ind(size(A),SI,mod(SJ+1,3)+1))];
    % Vertex opposite second subdivided edge
    [SJ,SI] = find(M'==1 & ~first);
    I2 = sub2ind(size(A),SI,SJ);
    E13 = [E13; ...
      F(sub2ind(size(A),SI,mod(SJ,3)+1)) ...
      F(sub2ind(size(A),SI,mod(SJ+1,3)+1))];

    % New vertex positions at midpoints
    U13 = (V(E13(:,1),:)+V(E13(:,2),:))/2;
    % indices of midpoints
    nu = size(U13,1);
    i1 = offset+(1:(nu/2))';
    i2 = offset+((nu/2)) + (1:(nu/2))';
    temp1 = i1;
    i1(flip) = i2(flip);
    i2(flip) = temp1(flip);
    F13 = [i1 O1 i2;i2 O2 O3;i2 O3 i1];
    %reshape(F13,[],9)
  end

  %      o           o     
  %     / \         /|\     
  %    /   \  ---> / | \    
  %   /     \     /  |  \  
  %  o---x---o   o---o---o  
  function [U12,F12,E12] = one_two(offset,V,F,M)
    [SJ,SI] = find(M'==1);
    O1  = F(sub2ind(size(M),SI,SJ));
    O2  = F(sub2ind(size(M),SI,mod(SJ+0,3)+1));
    O3  = F(sub2ind(size(M),SI,mod(SJ+1,3)+1));
    E12 = [O2 O3];
    % New vertex positions at midpoints
    U12 = (V(E12(:,1),:)+V(E12(:,2),:))/2;
    nu = size(U12,1);
    i1 = offset + (1:nu)';
    F12 = [O1 O2 i1;O1 i1 O3];
  end

  keep_duplicates = false;
  sel = [];
  iters = 1;
   % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'KeepDuplicates','OnlySelected','Iterations'}, ...
    {'keep_duplicates','sel','iters'});
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
  if isempty(sel)
    sel = (1:size(F,1))';
  end

  if iters<1
    FF = F;
    VV = V;
    FO = (1:size(F,1))';
    return;
  end

  if islogical(sel)
    sel = find(sel);
  end
  sel = sel(:);

  switch size(F,2)
  % http://mathoverflow.net/questions/28615/tetrahedron-splitting-subdivision
  case 3
    % Add a new vertex at the midpoint of each edge
    nsel = setdiff((1:size(F,1))',sel);

    Fsel = F(sel,:);

    E = [F(:,2) F(:,3);F(:,3) F(:,1);F(:,1) F(:,2)];
    Esel = [F(sel,2) F(sel,3);F(sel,3) F(sel,1);F(sel,1) F(sel,2)];
    Ensel = [F(nsel,2) F(nsel,3);F(nsel,3) F(nsel,1);F(nsel,1) F(nsel,2)];
    
    [I] = ismember(sort(Ensel,2),sort(Esel,2),'rows');
    nn = numel(nsel);
    M = sparse(repmat(1:nn',1,3),repmat([1 2 3],nn,1),reshape(I,nn,3),nn,3);
    C = sum(M,2);
    M13 = M(C==2,:);
    S13 = nsel(C==2);
    M12 = M(C==1,:);
    S12 = nsel(C==1);

    n = size(V,1);
    S14 = union(sel,nsel(C==3));

    S11 = setdiff((1:size(F,1))',[S14(:);S13(:);S12(:)]);
    [U14,F14,EU14] = one_four(size(V,1),V,F(S14,:));
    [U13,F13,EU13] = one_three(size(V,1)+size(U14,1),V,F(S13,:),M13);
    [U12,F12,EU12] = one_two(size(V,1)+size(U14,1)+size(U13,1),V,F(S12,:),M12);
    F11 = F(S11,:);

    FF = [F14;F13;F12;F11];
    FO = [S14;S14;S14;S14;S13;S13;S13;S12;S12;S11];
    U = [U14;U13;U12];
    EU = [EU14;EU13;EU12];
    nu = size(U,1);

    % find unique midpoints (potentially slow, though matlab is good at
    % these)
    if keep_duplicates
      U = U;
      J = (1:nu)';
    else
      [~,I,J] = unique(sort(EU,2),'rows');
      U = U(I,:);
    end
    % append unique midpoints to vertex positions
    VV = [V ; U];
    % reindex map from duplicate midpoint indices to unique midpoint indices
    J = [(1:n)';J+n];
    % reindex faces
    FF = J(FF);
  case 2
    if numel(sel)==size(F,1)
      m = [ (V(F(:,1),:) + V(F(:,2),:))/2 ];
      % indices of new midpoint vertices
      im = size(V,1) + (1:size(m,1))';
      % insert new face indices
      FF = [F(:,1) im;im F(:,2)];
      nf = size(F,1);
      FO = [1:nf 1:nf]';
      % append unique midpoints to vertex positions
      VV = [V;m];
      % No duplicates in 2D case
    else
      Fsel = F(sel,:);
      nsel = setdiff(1:size(F,1),sel);
      Fnsel = F(nsel,:);
      [VV,FF,FO] = upsample(V,Fsel,'KeepDuplicates',keep_duplicates);
      FF = [Fnsel;FF];
      FO = [nsel;sel(FO)];
    end
  end

  % Recursive call (iters=0 base case will be handled at top)
  if isempty(sel)
    sel = 1:size(F,1);
  end
  [VV,FF,FOr] = upsample( ...
    VV,FF, ...
    'OnlySelected',find(ismember(FO,sel)), ...
    'Iterations',iters-1, ...
    'KeepDuplicates',keep_duplicates);
  FO = FO(FOr);

end
