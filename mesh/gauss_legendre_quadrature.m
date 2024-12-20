function [y,u,x,w] = gauss_legendre_quadrature(N0,a,b)
  % [y,u,x,w] = gauss_legendre_quadrature(N0,a,b)
  % 
  % Inputs:
  %   nq  number of quadrature points
  %   a   lower bound of the integration
  %   b   upper bound of the integration
  % Outputs:
  %   y  Legendre polynomial zeros so that x = (a*(1-y)+b*(1+y))/2
  %   u  raw weights so that w = (b-a).*uw
  %   x  quadrature points
  %   w  quadrature weights
  %

  % Adapated from:
  % lgwt.m
  %
  % This script is for computing definite integrals using Legendre-Gauss 
  % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
  % [a,b] with truncation order N
  %
  % Suppose you have a continuous function f(x) which is defined on [a,b]
  % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
  % the values contained in the x vector to obtain a vector f. Then compute
  % the definite integral using sum(f.*w);
  %
  % Written by Greg von Winckel - 02/25/2004
  % slighly modified by Seungbae Bang - 01/18/2021
  % further modified by Alec Jacobson - 2024

  N=N0-1;
  N1=N+1; N2=N+2;
  
  xu=linspace(-1,1,N1)';
  
  % Initial guess
  y=-(cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2));
  
  % Legendre-Gauss Vandermonde Matrix
  L=zeros(N1,N2);
  
  % Derivative of LGVM
  Lp=zeros(N1,N2);
  
  % Compute the zeros of the N+1 Legendre Polynomial
  % using the recursion relation and the Newton-Raphson method
  
  y0=zeros(N+1,1)+2;
  
  % Iterate until new points are uniformly within epsilon of old points
  while max(abs(y-y0))>eps 
    L(:,1)=1;
    L(:,2)=y;
    for k=2:N1
      L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0=y;
    y=y0-L(:,N2)./Lp;
  end
  
  % Compute the weights
  u = (N2/N1)^2./((1-y.^2).*Lp.^2);

  if nargout > 2
    if nargin < 3
      b = 1;
      if nargin < 2
        a = 0;
      end
    end

    % Linear map from [-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;    
    w=(b-a).*u;
  end


end
