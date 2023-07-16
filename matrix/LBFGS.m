classdef LBFGS
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
      if nargin>=1
        this.m = m;
      end
    end

    % Alogorithm 7.5 Nocedal and Wright
    function [this,dx] = step(this,xk,dfdxk)
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

