function [phi,psi,N,UN] = green_coordinates(C,E,eta)
  % GREEN_COORDINATES  Compute green coordinates at evaluation points eta
  % w.r.t. cage (C,E)
  %
  % [phi,psi,N,UN] = green_coordinates(C,E,eta)
  %
  % Inputs:
  %   C  #C by dim list of cage vertices
  %   E  #E by dim list of cage facet indices into C
  %   eta  #eta by dim list of evaluation points
  % Outputs:
  %   phi  #eta by #C list of phi coordinates
  %   psi  #eta by #E list of psi coordinates
  %   N   #E by dim list of unit facet normals
  %
  % Example:
  %   % Cage in (C,E) and mesh in (V,F) _inside_ of (C,E)
  %   [phi,psi,N,UN] = green_coordinates(C,E,V);
  %   % Reproduce rest pose
  %   U = phi*C + psi*N;
  %   assert(max(abs(V-U))<eps);
  %

  % Following algorithm 1 in "Green Coordinates" [Lipman et al. 2008]
  A = C(E(:,2),:)-C(E(:,1),:);
  % B(i,:,k) = C(E(i,1),:) - V(k,:)
  B = bsxfun(@minus,C(E(:,1),:),permute(eta,[3 2 1]));
  Q = sum(A.*A,2);
  S = sum(B.*B,2);
  R = 2*sum(bsxfun(@times,A,B),2);
  % normal
  UN = A*[0 -1;1 0];
  N = normalizerow(UN);
  Anorm = normrow(A);
  BA = sum(bsxfun(@times,B,bsxfun(@times,Anorm,N)),2);
  X = bsxfun(@minus,4*bsxfun(@times,S,Q),R.^2);;
  X = max(X,0);
  SRT = sqrt(X);
  L0 = log(S);
  L1 = log(bsxfun(@plus,bsxfun(@plus,S,Q),R));
  A0 = atan2(R,SRT)./SRT;
  A1 = atan2(bsxfun(@plus,2*Q,R),SRT)./SRT;
  A10 = A1-A0;
  L10 = L1-L0;
  PHI2 = permute( ...
    -bsxfun(@times,BA/(2*pi), ...
      bsxfun(@rdivide,L10,2*Q)-bsxfun(@times,A10,bsxfun(@rdivide,R,Q))), ...
      [3 1 2]);
  PHI1 = permute( ...
    bsxfun(@times,BA/(2*pi), ...
      bsxfun(@rdivide,L10,2*Q)-bsxfun(@times,A10,2+bsxfun(@rdivide,R,Q))), ...
      [3 1 2]);
  % points lying on edges
  for e = 1:size(E,1) 
    [t,sqr_d] = project_to_lines(eta,C(E(e,1),:),C(E(e,2),:));
    on_edge = ((abs(sqr_d) < 2e-7) & ((t > -1e-10) & (t < (1+1e-10))));
    PHI1(on_edge,e) = -0.5*(1-t(on_edge));
    PHI2(on_edge,e) = -0.5*(t(on_edge));
    on_line = ((abs(sqr_d) < 2e-7) & ((t < -1e-10) | (t > (1+1e-10))));
    PHI1(on_line,e) = 0;
    PHI2(on_line,e) = 0;
  end
  phiI = repmat(1:size(eta,1),size(E,1)*2,1)';
  phiJ = repmat(E(:),1,size(eta,1))';
  phiV = [PHI1 PHI2];
  phi = -full(sparse(phiI,phiJ,phiV,size(eta,1),size(C,1)));
  % point at vertices
  EV = C(E(:,2),:)-C(E(:,1),:);
  A = full(sparse(E(:,[1 1]),repmat([1 2],size(E,1),1), EV,size(C,1),3));
  B = full(sparse(E(:,[2 2]),repmat([1 2],size(E,1),1),-EV,size(C,1),3));
  alpha = atan2(B(:,2),B(:,1))-atan2(A(:,2),A(:,1));
  alpha(alpha<0) = alpha(alpha<0)+2*pi;
  sqr_d = pdist2(eta,C).^2;
  [I,J] = find(sqr_d< 2e-7);
  phi(sub2ind(size(phi),I,J)) = (pi-alpha(J))/(2*pi)+0.5;

%phiC = phiC - diag(diag(phiC)) + diag((pi-alpha)/(2*pi)+0.5);

  psi = ...
    bsxfun(@times,-Anorm/(4*pi), ...
      bsxfun(@times,bsxfun(@minus,4*S,bsxfun(@rdivide,R.^2,Q)),A10)+ ...
      bsxfun(@plus,bsxfun(@times,bsxfun(@rdivide,R,2*Q),L10),L1)-2);
  psi = permute(psi,[3 1 2]);
  % points lying on edges
  for e = 1:size(E,1) 
    [t,sqr_d] = project_to_lines(eta,C(E(e,1),:),C(E(e,2),:));
    %on_edge = ((abs(sqr_d) < 2e-7) & ((t > -1e-10) & (t < (1+1e-10))));
    on_edge = abs(sqr_d) < 2e-7;
    % edge length
    s = norm(C(E(e,2),:)-C(E(e,1),:));
    % From maple
    toe = s*t(on_edge);
    psi(on_edge,e) = ((toe-s).*log((toe-s).^2)-toe.*log(toe.^2)+2*s)/(4*pi);
    on_end = on_edge & (abs(t)<1e-8 | abs(1-t)<1e-8);
    psi(on_end,e) = (-s*log(s^2)+2*s)/(4*pi);
  end

end
