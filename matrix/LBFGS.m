classdef LBFGS
  % LBFGS  Class implementing an iteration step direction for the LBFGS algorithm.
  % This class prepares the step direction dx and the caller is responsible for
  % conducting the step (e.g., using backtracking_line_search).
  %
  % L = LBFGS() creates an LBFGS object with a default number {10} of past states.
  % L = LBFGS(m) creates an LBFGS object with a memory of m past states.
  %
  % LBFGS methods:
  %    step - compute the step direction
  % 
  % Example:
  %   f =    @(x) 0.26*(x(1)^2+x(2)^2) - 0.48*x(1)*x(2);
  %   grad = @(x) [0.52*x(1)-0.48*x(2); 0.52*x(2)-0.48*x(1)];
  %   L = LBFGS(2);
  %   x = [0;10];
  %   for iter = 1:10
  %     dfdx = grad(x);
  %     % compare to `dx = -dfdx;`
  %     [L,dx] = step(L,x,dfdx);                                  
  %     [t,x] = backtracking_line_search(f,x,dfdx,dx,0.001,0.5,100);
  %     fprintf(['%2d: f = %0.1g\n'],iter,f(x));
  %   end

  properties
    % number of past states to remember
    m = 10
    % current interation
    k = 0
    % past state
    xk = []
    dfdxk = []
    s = []
    y = []
    rho = []
    alpha = []
  end
  methods

    function this = LBFGS(m)
      % this = LBFGS(m)
      if nargin>=1
        this.m = m;
      end
    end

    % Alogorithm 7.5 Nocedal and Wright
    function [this,dx] = step(this,xk,dfdxk)
      % LBFGS.STEP Compute the step direction given the stored state and new
      % value and gradients.
      %
      % L = LBFGS();
      % …
      % [L,dx] = step(L,xk,dfdxk)
      %
      % Inputs:
      %   xk  #x array of current values for iteration k
      %   dfdxk  #x array of current gradients for iteration k
      % Outputs:
      %   dx  #x array of step direction values
      % 

      % Normalize inputs to vectors
      xk_size = size(xk);
      xk = reshape(xk,[],1);
      dfdxk = reshape(dfdxk,[],1);

      % convert iteration to limited memory index
      iter2ind = @(i) mod(i-1,this.m)+1;

      % Increment state
      if this.k == 0
        n = numel(xk);
        this.s = nan(n,this.m);
        this.y = nan(n,this.m);
        this.rho = nan(1,this.m);
        this.alpha = nan(1,this.m);
      else 
        K = iter2ind(this.k);
        this.s(:,K) =    xk - this.xk;
        this.y(:,K) = dfdxk - this.dfdxk;
        this.rho(K) = 1/(this.y(:,K)'*this.s(:,K));
      end
      this.xk = xk;
      this.dfdxk = dfdxk;
      this.k = this.k + 1;
      % we are now computing dx = -Hₖ ∇fₖ
      if this.k==1
        dx = reshape(-dfdxk,xk_size);
        return
      end

      % Alogorithm 7.4 Nocedal and Wright
      q = dfdxk;
      for i = this.k-1:-1:max(this.k-this.m,1)
        I = iter2ind(i);
        this.alpha(I) = this.rho(I)*(this.s(:,I)'*q);
        q = q - this.alpha(I)*this.y(:,I);
      end
      if nargin<4
        K_minus_1 = iter2ind(this.k-1);
        gamma_k =  ...
          (this.s(:,K_minus_1)'*this.y(:,K_minus_1)) / ...
          (this.y(:,K_minus_1)'*this.y(:,K_minus_1));
        H0k = gamma_k; % times identity; just keep scalar.
      end
      r = H0k * q;
      for i = max(this.k-this.m,1):this.k-1
        I = iter2ind(i);
        beta = this.rho(I)*(this.y(:,I)'*r);
        r = r + this.s(:,I)*(this.alpha(I)-beta);
      end
      dx = reshape(-r,xk_size);
    end

  end
end

