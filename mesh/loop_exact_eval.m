function [x,S,B0,k] = loop_exact_eval(V,F,f,bc,AF,varargin)
  % LOOP_EXACT_EVAL Exactly evaluate the limit position of point on a a Loop
  % subdivision surface (V,F) specified by the face f and barycentric
  % coordinates bc
  %
  %[x,S,B0,k] = loop_exact_eval(V,F,f,bc)
  %[x,S,B0,k] = loop_exact_eval(V,F,f,bc,AF)
  %
  % Inputs:
  %   V  #V by dim list of vertex coordinates
  %   F  #F by 3 list of triangles indices into V
  %   f  an index in F of face containing query 
  %   bc  3-vector of barycentric coordinates of query in face f
  %   AF  #F by #F sparse adjacency matrix, so that AF(i,j) is nonzero unless
  %     face i and j share a **vertex** {[]}
  %%   Optional:
  %%     'FaceFaceAdjacency' followed by AF so that AF(i,j) is nonzero unless
  %%       face i and j share a **vertex** {[]}.
  % Outputs:
  %   x  dim-long limit position of equery
  %   S  12 by #V sparse matrix containing subdivision coeffients collected
  %      down to regular patch containing query. The 12 rows are in Stam's
  %      order
  %   B0  12-long vector of basis functions at regular patch evaluated at query
  %     (same order as S)
  %   k  number of subdivisions needed to reach regular patch
  %
  % see also: stam_order
  %

  %AF = [];
  %params_to_variables = containers.Map( ...
  %  {'FaceFaceAdjacency'}, ...
  %  {'AF'});
  %v = 1;
  %while v <= numel(varargin)
  %  param_name = varargin{v};
  %  if isKey(params_to_variables,param_name)
  %    assert(v+1<=numel(varargin));
  %    v = v+1;
  %    % Trick: use feval on anonymous function to use assignin to this workspace
  %    feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
  %  else
  %    error('Unsupported parameter: %s',varargin{v});
  %  end
  %  v=v+1;
  %end

  % Well, actually, interior edges are just fine. Interior regular vertices are
  % also just. Boundary edges should be fine, but we'd have to implement the
  % reduction to 1D splines. Then boundary vertices would also be just fine. So
  % **_really_** this _should_ just be an assert on queries at irregular
  % vertices.
  %assert(~any(bc<=0),'query cannot lie on/past a vertex or edge');
  assert(~any(bc>=1),'query cannot lie on/past vertex');
  if nargin < 5
    AF = [];
  end
  bc0 = bc;
  f0 = f;
  V0 = V;
  F0 = F;
  S = speye(size(V,1));
  k = 0;
  while true
    %tsurf(F0,V0,'FaceIndices',1,'CData',sparse(f0,1,1,size(F0,1),1))
    %set(gca,'YDir','reverse');
    %pause

    % limit F to f and f's vertex-neighbors
    if isempty(AF)
      FV = sparse(repmat(1:size(F0,1),3,1)',F0,1,size(F0,1),size(V0,1));
      AF = FV*FV';
    end
    N = unique([f0;find(AF(:,f0))]);
    f0 = find(f0==N);
    F0 = F0(N,:);
    [V0,UI,UJ] = remove_unreferenced(V0,F0);
    S = S(UJ,:);
    % reshape for stupid single triangle case
    F0 = reshape(UI(F0),size(F0));
    % Now, (V0,F0) is a tight patch around f0

    % For each edge touching the query, check if it's a boundary edge. These
    % are not handled, though they could be. We only need to check this for the
    % coarsest level because to appear on the boundary later, the query must
    % have been on the boundary to begin with.
    if k==0
      bc_edges = find(bc<=0);
      if ~isempty(bc_edges)
        [~,B] = on_boundary(F0);
        Bf0 = B(f0,:);
        assert(~any(Bf0(bc_edges)),'Query cannot lie on/past boundary edge');
      end
    end

    %tsurf(F0,V0,'FaceIndices',1,'CData',sparse(f0,1,1,size(F0,1),1))
    %set(gca,'YDir','reverse');
    %pause

    % determine that mesh is regular around f, O(N)
    A = adjacency_matrix(F0);
    % V0 should include 12 vertices and F0 should include 13 faces and all
    % vertices of f have valence 6
    is_reg = ...
      size(V0,1) == 12 && size(F0,1) == 13 && all(sum(A(:,F0(f0,:)),1) == 6);

    if is_reg
      % determine reordering so that face f is the [4 7 8] face of stam's
      % regular patch.
      I = stam_order(F0,A);
      % Exact evaluation of triangle spline functions for Jos Stam's regular
      % patch: vertex order matters.
      Bfun = @(v,w) [ ...
        v.*((v+w-1).^(3)).*(-1/6)+((v+w-1).^(4)).*(1/12), ...
        w.*((v+w-1).^(3)).*(-1/6)+((v+w-1).^(4)).*(1/12), ...
        (v.*v.*v).*(v.*6+w.*6-6).*(-1/12)-v.*((v+w-1).^(3)).*(1/2)-w.*((v+w-1).^(3)).*(1/6)+(v.*v.*v).*w.*(1/6)+((v+w-1).^(4)).*(1/12)+(v.*v.*v.*v).*(1/12)+(v.*v).*((v+w-1).^(2))-(v.*v).*w.*(v.*6+w.*6-6).*(1/12)+v.*w.*((v+w-1).^(2)).*(1/2), ...
        (v.*v).*(w.*w)-(v.*v.*v).*(v.*8+w.*8-8).*(1/12)-(w.*w.*w).*(v.*8+w.*8-8).*(1/12)-v.*((v+w-1).^(3)).*2-w.*((v+w-1).^(3)).*2+v.*(w.*w.*w).*(1/2)+(v.*v.*v).*w.*(1/2)+((v+w-1).^(4)).*(1/2)+(v.*v.*v.*v).*(1/12)+(w.*w.*w.*w).*(1/12)+(v.*v).*((v+w-1).^(2)).*2+(w.*w).*((v+w-1).^(2)).*2-v.*(w.*w).*(v.*36+w.*36-36).*(1/12)-(v.*v).*w.*(v.*36+w.*36-36).*(1/12)+v.*w.*((v+w-1).^(2)).*5, ...
        (w.*w.*w).*(v.*6+w.*6-6).*(-1/12)-v.*((v+w-1).^(3)).*(1/6)-w.*((v+w-1).^(3)).*(1/2)+v.*(w.*w.*w).*(1/6)+((v+w-1).^(4)).*(1/12)+(w.*w.*w.*w).*(1/12)+(w.*w).*((v+w-1).^(2))-v.*(w.*w).*(v.*6+w.*6-6).*(1/12)+v.*w.*((v+w-1).^(2)).*(1/2), ...
        (v.*v.*v).*(v.*2+w.*2-2).*(-1/12)+(v.*v.*v.*v).*(1/12), ...
        (v.*v).*(w.*w).*2-(v.*v.*v).*(v.*24+w.*24-24).*(1/12)-(w.*w.*w).*(v.*6+w.*6-6).*(1/12)-v.*((v+w-1).^(3)).*(2/3)-w.*((v+w-1).^(3)).*(1/2)+v.*(w.*w.*w).*(2/3)+(v.*v.*v).*w.*2+((v+w-1).^(4)).*(1/12)+(v.*v.*v.*v).*(1/2)+(w.*w.*w.*w).*(1/12)+(v.*v).*((v+w-1).^(2)).*2+(w.*w).*((v+w-1).^(2))-v.*(w.*w).*(v.*36+w.*36-36).*(1/12)-(v.*v).*w.*(v.*6E1+w.*6E1-6E1).*(1/12)+v.*w.*((v+w-1).^(2)).*3, ...
        (v.*v).*(w.*w).*2-(v.*v.*v).*(v.*6+w.*6-6).*(1/12)-(w.*w.*w).*(v.*24+w.*24-24).*(1/12)-v.*((v+w-1).^(3)).*(1/2)-w.*((v+w-1).^(3)).*(2/3)+v.*(w.*w.*w).*2+(v.*v.*v).*w.*(2/3)+((v+w-1).^(4)).*(1/12)+(v.*v.*v.*v).*(1/12)+(w.*w.*w.*w).*(1/2)+(v.*v).*((v+w-1).^(2))+(w.*w).*((v+w-1).^(2)).*2-(v.*v).*w.*(v.*36+w.*36-36).*(1/12)-v.*(w.*w).*(v.*6E1+w.*6E1-6E1).*(1/12)+v.*w.*((v+w-1).^(2)).*3, ...
        (w.*w.*w).*(v.*2+w.*2-2).*(-1/12)+(w.*w.*w.*w).*(1/12), ...
        (v.*v.*v).*w.*(1/6)+(v.*v.*v.*v).*(1/12), ...
        (v.*v).*(w.*w)-(v.*v.*v).*(v.*2+w.*2-2).*(1/12)-(w.*w.*w).*(v.*2+w.*2-2).*(1/12)+v.*(w.*w.*w).*(1/2)+(v.*v.*v).*w.*(1/2)+(v.*v.*v.*v).*(1/12)+(w.*w.*w.*w).*(1/12)-v.*(w.*w).*(v.*6+w.*6-6).*(1/12)-(v.*v).*w.*(v.*6+w.*6-6).*(1/12), ...
        v.*(w.*w.*w).*(1/6)+(w.*w.*w.*w).*(1/12), ...
      ];
      B0 = Bfun(bc0(2),bc0(3));
      S0 = speye(12);
      % reorder from stams ordering
      S0 = S0(I,:);
      % In the regular case:
      S = S0*S;
      break;
    else
      [V1,F1,S1,SJ] = loop(V0,F0,1);
      S = S1*S;
      % which child?
      ch = (bc0(1)>1/2)*1 + (bc0(2)>1/2)*2 + (bc0(3)>1/2)*3 + all(bc0<=1/2)*4;
      f1 = f0+(ch-1)*size(F0,1);
      % New barys
      bc1 = (ch==4)+(-(ch==4)*2+1)*bc0*2-(ch==(1:3));
      % Roll
      bc1 = bc1([mod(ch-1,3)+1:3 1:mod(ch-1,3)]);
      % get ready for recursion
      bc0 = bc1;
      F0 = F1;
      V0 = V1;
      f0 = f1;
      k = k+1;
      AF = [];
    end
  end
  x = B0 * S * V;
end
