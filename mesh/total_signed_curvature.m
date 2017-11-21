function [tsc,alpha,V] = total_signed_curvature(P,E)
  % TOTAL_SIGNED_CURVATURE  Total discrete signed curvature of a given closed
  % planar curve P (see http://ddg.cs.columbia.edu/SIGGRAPH05/Didactic.pdf)
  %
  % [tsc,alpha] = total_signed_curvature(P)
  % [tsc,alpha] = total_signed_curvature(P,E)
  %
  % Inputs:
  %   P  #P by 2 list of positions of a closed piecewise-linear curve
  %   or 
  %   P  #P by 2 list of positions 
  %   E  #E by 2 list of orientedd edge indices into P
  % Outputs:
  %   tsc  total signed curvature scalar, should be 2*pi*k where k∈0,1,2,...
  %   alpha  #P by 1 list of exterior angles (boundaries will get 0)
  %   V  #E by 2 list of edge vectors
  %

  if nargin==1
    E = [1:size(P,1);2:size(P,1) 1]';
  end

  % edge vectors
  V = P(E(:,2),:)-P(E(:,1),:);

  % Incoming and outgoing vectors (this will produce non-sense on boundaries and
  % non-manifold curves
  IV = sparse(E(:,[1 1]),repmat([1 2],size(E,1),1),V,size(P,1),2);
  OV = sparse(E(:,[2 2]),repmat([1 2],size(E,1),1),V,size(P,1),2);


  % Exterior angles between consecutive edge vectors
  alpha = atan2( ...
    IV(:,1).*OV(:,2) - IV(:,2).*OV(:,1), ...
    IV(:,1).*OV(:,1) + IV(:,2).*OV(:,2) );

  alpha(isnan(alpha) | isinf(alpha)) = 0;
  
  % total signed curvature
  tsc = sum(alpha);

end
