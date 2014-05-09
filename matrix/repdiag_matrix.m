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

    function C = mtimes(this,that)
      % MTIMES Computes C = this * that where this and/or that is a repeated
      % block diagonal matrix.
      %
      % Inputs:
      %   this  n by m matrix of some form (repdiag_matrix or double)
      %   that  m by p matrix of some form (repdiag_matrix or double)
      % Outputs:
      %   C  n by p matrix
      %
      % Compare to:
      %   C = this.expand() * that.expand();
      %   C = this * that.expand();
      %   C = this.expand() * that;

      if isa(this,'repdiag_matrix') && ~isa(that,'repdiag_matrix')
        m = size(this.A,1);
        B = that;
        assert(size(this.A,2)*this.d == size(B,1));
        n = size(B,1)/this.d;
        assert(n == m);
        k = size(B,2);
        % Apply to columns of each block (disregarding order but then
        % reassembling)
        BB = reshape(B,n,k*this.d);
        CC = this.A*BB;
        C = reshape(CC,m*this.d,k);
      elseif ~isa(this,'repdiag_matrix') && isa(that,'repdiag_matrix')
        %% Cheapskate
        %C = (that'*this')';
        B = this;
        % B * this
        % k by n*d  *  m*d by m*d
        m = size(that.A,1);
        assert(size(that.A,2)*that.d == size(B,2));
        n = size(B,2)/that.d;
        assert(n == m);
        k = size(B,1);
        % Apply to columns of each block (disregarding order but then
        % reassembling)
        % STILL NEED TO TRANSPOSE... Lame...
        BB = reshape(B',n,k*that.d);
        % BB * A
        % k*d by n  *  n * n
        CC = that.A'*BB;
        C = reshape(CC,m*that.d,k)';
      else
        assert(this.d == that.d, ...
          'Non-matching repdiag_matrix multiplication not supported');
        assert(all(size(this.A)==size(that.A)), ...
          'Non-matching repdiag_matrix multiplication not supported');
        C = repdiag_matrix(this.A*that.A,this.d);
      end

    end

    function C = mldivide(this,B)
      % MLDIVIDE Computes C = this \ B where this is a repeated block diagonal
      % matrix and B is a regular old matrix
      %
      % Inputs:
      %   this  d*m by d*n repeated block diagonal matrix
      %   B  d*n by k matrix
      % Outputs:
      %   C  d*m by k matrix
      %
      % Compare to:
      %   C = A.expand() \ B;

      assert(size(this.A,2)*this.d == size(B,1));
      n = size(B,1)/this.d;
      k = size(B,2);
      m = size(this.A,1);
      
      %C = this.expand() * B;
      BB = reshape(B,n,k*this.d);
      CC = this.A \ BB;
      C = reshape(CC,m*this.d,k);
    end

    function varargout = chol(varargin)
      % CHOL compute cholesky decompostion of repeated block matrix
      %
      % See also: chol
      this = varargin{1};
      varargin{1} = this.A;
      switch nargout
      case 1
        [varargout{1}] = chol(varargin{:});
      case 2
        [varargout{1},varargout{2}] = chol(varargin{:});
      case 3
        [varargout{1},varargout{2},varargout{3}] = chol(varargin{:});
      otherwise
        error(['Too many outputs (' num2str(nargout) ')']);
      end
      varargout{1} = repdiag_matrix(varargout{1},this.d);
      varargout{3} = repdiag_matrix(varargout{3},this.d);
    end

    function T = transpose(this)
      % TRANSPOSE computes transpose of block diagonal matrix
      T = repdiag_matrix(transpose(this.A),this.d);
    end

    function T = ctranspose(this)
      % CTRANSPOSE computes conjugate transpose of block diagonal matrix
      T = repdiag_matrix(ctranspose(this.A),this.d);
    end

  end
  
end
