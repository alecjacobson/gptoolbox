function [Q,F] = resample_polygon(P,E,h)
  % RESAMPLE_POLYGON  Sample edges of a polygon so all new edges are smaller
  % than h
  %
  % Inputs:
  %   P  #P by dim list of input points
  %   E  #E by 2 list of input edge indices into P
  %   h  upper bound on edge length
  % Outputs:
  %   Q  #Q by dim list of output points
  %   F  #F by 2 list of output edge indices into Q
  %

  Q = P;
  F = [];
  l = edge_lengths(P,E);
  for e = 1:size(E,1)
    if l(e)<=h
      % edge is already short enough
      F = [F;E(e,:)];
      continue;
    end
    % source and dest
    s = E(e,1);
    Ps = P(s,:);
    d = E(e,2);
    Pd = P(d,:);
    num_segs = ceil(l(e)/h);
    t = linspace(0,1,num_segs+1)';
    M = bsxfun(@plus,Ps,bsxfun(@times,t(2:end-1),(Pd-Ps)));
    F = [F; [s size(Q,1)+(1:size(M,1)); size(Q,1)+(1:size(M,1)) d]'];
    Q = [Q;M];
  end
end
