function [kappa,alpha,ev,l] = curvature(P,varargin)
  % CURVATURE compute pointwise curvature on a closed piecewise-linear curve P
  % 
  % [kappa,alpha,ev,l] = curvature(P)
  % [kappa,alpha,ev,l] = curvature(P,E)
  %
  % Inputs:
  %   P  #P by 2 list of curve vertices
  %    or
  %   P  #P by 2 list of curve vertices
  %   E  #E by 2 list of directed curve edges into P
  % Outputs:
  %   kappa  #P list of curvature values
  %   alpha  #P list of exterior angle values
  %   ev  #P by 2 list of edge vectors, ev(i,:) goes from P(i,:) to P(i+1,:)
  %   l  #P list of edge vector lengths, ev(i,:) goes from P(i,:) to P(i+1,:)
  %

  n = size(P,1);

  switch nargin
  case 1
    E = [1:n;2:n 1]';
  case 2
    E = varargin{1};
  end
  % edge vectors, evi goes from Pi to Pi+1
  ev = P(E(:,2),:) - P(E(:,1),:);

  % find incoming and outgoing edge of each vertex
  A = sparse(E(:,2),E(:,1),1,n,n);
  A = A-A';
  [~,EI] = max(A,[],2);
  [~,EO] = min(A,[],2);
  evI =  P - P(EI,:);
  evO =  P(EO,:) - P;

  % Exterior angles between consecutive edge vectors
  alpha = atan2( ...
    evI(:,1).*evO(:,2) - evI(:,2).*evO(:,1), ...
    evI(:,1).*evO(:,1) + evI(:,2).*evO(:,2) );

  lI = sqrt(sum(evI.^2,2));
  lO = sqrt(sum(evO.^2,2));
  kappa = alpha ./ (0.5 * (lI + lO));
end
