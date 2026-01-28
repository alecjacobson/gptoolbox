function [E,PI,PC] = cells(d1,d2)
  % [E,PI,PC] = cells(d1,d2)
  %
  % Inputs:
  %   d1  #V by #E sparse, oriented boundary matrix
  %   d2  #E by #F sparse, oriented boundary matrix
  % Outputs:
  %   E  #E by 2 list of 1-cells indexing 1:#V, except blank edges will be [0 0]
  %   PI  #PI streamed list of polygon indices into rows of V
  %   PC  #PC+1 streamed list of cumulative counts of polygon indices
  %
  % see also: polygons_to_triangles, boundary_matrix
  %
  % Example:
  %   [E,PI,PC] = cells(d1,d2);
  %   F = polygons_to_triangles(PI,PC);
  %   % or if you know nat all(sum(abs(d2)>0)) == 3
  %   F = reshape(PI,3,[]).';


  % size of each 2-cell
  K = full(sum(abs(d2)>0)).';
  blank_edges = all(d1==0,1);
  assert(all(sum(d1(:,~blank_edges)==-1) == 1),'every edge has exactly one -1');
  assert(all(sum(d1(:,~blank_edges)==+1) == 1),'every edge has exactly one +1');
  [valid_I,EI] = max(d1<0);
  [valid_J,EJ] = max(d1>0);
  valid = valid_I & valid_J;
  assert(all(valid| (~valid_I & ~valid_J)));
  E = [EI;EJ]';
  E(~valid,:) = 0;
  if nargout == 1
    return;
  end

  % half-edge
  h1 = [d1 -d1];
  % half-edge 2-cell incidence
  h2 = [d2>0;d2<0];

  %S = [speye(size(d2,1)),-speye(size(d2,1))];
  %% d2 = S*h2

  % hn(i) = j indicates j is "next" half-edge of i.
  % hn(i) = 0 means i is not on a 2-cell (boundary)
  %   (h1'*(h1==-1)) possible next candidates
  %    AND
  %   (h2*h2') sharing a face
  Q = (h1'*(h1==-1)).*(h2*h2');
  [nv,hn] = max(Q==1,[],2);

  % mask out boundaries
  hn = hn.*nv;

  % And this isn't handling boundaries yet...
  check_simply_connected = false;
  if check_simply_connected
    H = sparse((1:numel(hn))',hn,1,numel(hn),numel(hn));
    H = H + H';
    [nc,C] = conncomp(H);
    C2h = sparse(C,1:numel(hn),1,nc,numel(hn));
    % ncp(i) = number of components of boundary edges in face i
    ncp = sum(full((C2h*abs(h2))>0));
    assert(all(ncp==1),'Every face should be simply connected (one component of boundary edges)');
  end

  % pick one half-edge per face
  [~,hi] = max(h2);
  hi = hi';
  % origin vertex of each half-edge
  [~,hv] = max(h1==-1);
  hv = hv';

  PC = cumsum([0;K]);
  PI = zeros(PC(end),1);
  % Which faces are still unfinished.
  I = (1:size(d2,2))';
  c = 0;
  while true
    % There might be zero length cells
    keep = K(I) > c;
    I = I(keep);
    hi = hi(keep);
    c = c + 1;
    if isempty(I)
      break;
    end

    PI( PC(I) + c ) = hv(hi);
    hi = hn(hi);
  end

end
