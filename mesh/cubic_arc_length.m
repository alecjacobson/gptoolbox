function [l,dlda,dldb,dldC,d2lda2,d2ldb2,d2ldadb] = cubic_arc_length(C,nq,a,b)
  %
  % [l,dlda,dldb] = cubic_arc_length(C,nq)
  % [l,dlda,dldb] = cubic_arc_length(C,[x w])
  % [l,dlda,dldb] = cubic_arc_length(C,nq,a,b)
  % [l,dlda,dldb] = cubic_arc_length(C,[y u],a,b)
  %
  % Inputs:
  %   C  4 by dim by 1|#a list of cubic control points
  %   nq  number of quadrature points or quadrature points and weights
  %     or
  %   [x w]  nq by 2 list of quadrature points and weights
  %     or
  %   [y u]  nq by 2 list of Legendre polynomial zeros so that x = (a*(1-y)+b*(1+y))/2
  %     and raw weights so that w = (b-a).u*
  %   a  #a list of integral interval start values
  %   b  #a list of integral interval end values
  % Outputs:
  %   l  #a list of arc lengths
  %   dlda  #a list of derivatives of arc lengths wrt a
  %   dldb  #a list of derivatives of arc lengths wrt b
  %
  % See also: gauss_legendre_quadrature
  %   
  % Example:
  %   l = cubic_arc_length(permute(reshape(P(C',:),4,size(C,1),size(P,2)),[1 3 2]));
  if nargin<2
    nq = 10;
  end
  if nargin<3
    a = 0;
    b = 1;
  end

  assert(numel(a) == numel(b));
  if size(nq,2) == 1
    assert(numel(nq) == 1);
    if numel(a) == 1
      [y,u,x,w] = gauss_legendre_quadrature(nq,a,b);
    else
      [y,u] = gauss_legendre_quadrature(nq);
      x = [];
      w = [];
    end
  elseif size(nq,2) == 2
    if nargin<3
      % Evaluation points and weights given directly
      x = nq(:,1);
      w = nq(:,2);
    else
      % canonical evaluation points and raw weights
      y = nq(:,1);
      u = nq(:,2);
      x = [];
      w = [];
    end
    nq = size(y,1);
  end

  if isempty(x)
    assert(isempty(w));
    x = (a.'.*(1-y)+b.'.*(1+y))/2;
    w = (b-a).'.*u;
  end
  x = permute(x,[1 3 2]);
  w = permute(w,[1 3 2]);

  C21 = C(2,:,:)-C(1,:,:);
  C32 = C(3,:,:)-C(2,:,:);
  C43 = C(4,:,:)-C(3,:,:);
  b21 = 3*(1-x).^2;
  b32 = 6*(1-x).*x;
  b43 = 3*x.^2;
  T = b21.*C21 + ...
      b32.*C32 + ...
      b43.*C43;
  s = sqrt(sum(T.^2,2));
  l = sum(w.*s,1);
  l = permute(l,[3 1 2]);

  if nargout>1
    dwda = -u;
    dwdb =  u;
    % dlda = d/da( wᵀ s)  = dw/daᵀ s + wᵀ ds/da
    dsdT = 1./s.*T;
    x6 = 6*x;
    dTdx = x6.*(C43-C32) + C21.*(x6-6) - C32.*(x6-6);
    dxda = 1/2 - y/2;
    dxdb = 1/2 + y/2;
    dsdx = sum(dsdT.*dTdx,2);
    dsda = dsdx.*dxda;
    dsdb = dsdx.*dxdb;
    dlda = sum(dwda.*s,1) + sum(w.*dsda,1);
    dldb = sum(dwdb.*s,1) + sum(w.*dsdb,1);
    dlda = permute(dlda,[3 1 2]);
    dldb = permute(dldb,[3 1 2]);

    if nargout>3
      % dldC = d/dC (wᵀ s) = dw/dCᵀ s + wᵀ ds/dC  = wᵀ ds/dC
      dTdC = zeros([nq size(C,2) numel(a) size(C)]);
      dTdC(:,1,:,:,1) = cat(4,-b21,b21-b32,b32-b43,b43);
      dTdC(:,2,:,:,2) = dTdC(:,1,:,:,1);
      dsdC = ((1./s).*sum(T.*dTdC,2));
      dldC = permute(sum(w.*dsdC,[1 2]),[3 4 5 1 2]);
    end

    if nargout>4
      % d²l/da² = d/da( dw/daᵀ s + wᵀ ds/da )
      %          = d²w/da²ᵀ s + 2 dw/daᵀ ds/da + wᵀ d²s/da²
      %          =       0ᵀ s + 2 dw/daᵀ ds/da + wᵀ d²s/da²
      d2wda2 = zeros(numel(u),1,numel(a));
      d2wdb2 = zeros(numel(u),1,numel(a));
      % dsda = sum(dsdT.*dTdx,2).*dxda;
      % dsda = dsdx.*dxda;
      % d/da(dsda) = d/da(dsdx) dxda + dsdx d/da(dxda)
      % d/da(dsda) = d/da(dsdx) dxda + dsdx 0
      d2Tdxda =  3*(y - 1).*(C(1,:,:) - 3*C(2,:,:) + 3*C(3,:,:) - C(4,:,:));
      d2Tdxdb = -3*(y + 1).*(C(1,:,:) - 3*C(2,:,:) + 3*C(3,:,:) - C(4,:,:));

      % dsdx = sum(dsdT.*dTdx,2);
      % dsdT = 1./s.*dT;
      % d/da(dsdT) = d/da(1./s.*dT) = -1./s² ds/da dT + 1./s d/da(dT)
      d2sdTda = -1./s.^2.*dsda.*T + 1./s.*dTdx.*dxda;
      d2sdTdb = -1./s.^2.*dsdb.*T + 1./s.*dTdx.*dxdb;

      d2sdxda = sum(d2sdTda.*dTdx,2) + sum(dsdT.*d2Tdxda,2);
      d2sdxdb = sum(d2sdTdb.*dTdx,2) + sum(dsdT.*d2Tdxdb,2);

      d2sda2 = d2sdxda.*dxda;
      d2sdb2 = d2sdxdb.*dxdb;
      d2lda2 = 2*sum(dwda.*dsda,1) + sum(w.*d2sda2,1);
      d2ldb2 = 2*sum(dwdb.*dsdb,1) + sum(w.*d2sdb2,1);

      % d²l/dadb = d/da( dw/dbᵀ s + wᵀ ds/db )
      %          = d²w/dadbᵀ s + dw/dbᵀ ds/da + dw/daᵀ ds/db + wᵀ d²s/dadb
      %          =       0ᵀ s + dw/dbᵀ ds/da + dw/daᵀ ds/db + wᵀ d²s/dadb
      % dsda = dsdx.*dxda;
      % d/db(dsda) = d/db(dsdx) dxda + dsdx d/db(dxda)
      % d/db(dxda) = 0
      % d/db(dsda) = d/db(dsdx) dxda + dsdx 0
      %            = d/db(dsdx) dxda
      d2sdadb = d2sdxdb .* dxda;
      d2ldadb = sum(dwdb.*dsda,1) + sum(dwda.*dsdb,1) + sum(w.*d2sdadb,1);
      d2lda2 = permute(d2lda2,[3 1 2]);
      d2ldb2 = permute(d2ldb2,[3 1 2]);
      d2ldadb = permute(d2ldadb,[3 1 2]);
    end

  end


end
