classdef repdiag_matrix
  % REPDIAG_MATRIX class for storing and dealing with matrices of the type:
  % M = repdiag(A);
  % aka M = diag(A,A,...,A);
  % aka M = kron(eye(d),A)
  % aka M is a block diagonal matrix with d copies of A along the diagonal
  %
  % Inputs:
  %   A  m by n matrix repeated along diagonal
  %   d  number of copies of A along diagonal
  % Output:
  %   this  an object giving access to the repeated matrix and fields making
  %     typical matrix operations faster by taking advantage of repeated block
  %     diagonal structure
  %
  % See also: repdiag

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read only class fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties(GetAccess=public,SetAccess=protected)
    % the m by n matrix to be repeated
    A 
    % the number of times A is repeated
    d 
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Private class fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties(Access=private)
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Public Class methods
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  methods(Access=public)

    function this = repdiag_matrix(varargin)
    % See also: repdiag_matrix, repdiag
      this.A = varargin{1};
      this.d = varargin{2};
    end

    function B = expand(this)
      % EXPAND Computes the true matrix represented by this class (i.e. calls
      % and returns repdiag(A,d)
      %
      B = repdiag(this.A,this.d);
    end

    function C = mtimes(this,B)
      % MTIMES Computes C = this * B where this is a repeated block diagonal
      % matrix and B is a regular old matrix
      %
      % Inputs:
      %   this  d*m by d*n repeated block diagonal matrix
      %   B  d*n by k matrix
      % Outputs:
      %   C  d*m by k matrix
      %
      % Compare to:
      %   C = A.expand() * B;
      
      %% naive way
      %C = this.expand() * B;

      % try this:
      % http://blogs.mathworks.com/videos/2012/01/10/how-to-do-a-matrix-reshape-by-blocks/
      assert(size(this.A,2)*this.d == size(B,1));
      n = size(B,1)/this.d;
      k = size(B,2);
      m = size(this.A,1);
      %C = cell(this.d,1);
      %for ii = 0:(this.d-1)
      %  C{ii+1} = this.A * B(ii*n + (1:n),:);
      %end
      %C = cell2mat(C);
      % Transpose each "block" of B, so that if B=[B1;B2;...;Bd];
      % now B = [B1 B2 ... Bd];
      BB = reshape(permute(reshape(B',[k n this.d]),[1 3 2]),[this.d*k n]);
      CC = BB * this.A';
      % now we have made C = [A*B1 A*B2 ... A*Bd]
      C = reshape(permute(reshape(CC,[k this.d n]),[1 3 2]),[k n*this.d])';

    end
  end
  
end
