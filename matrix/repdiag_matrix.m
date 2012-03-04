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

      assert(size(this.A,2)*this.d == size(B,1));
      n = size(B,1)/this.d;
      k = size(B,2);
      m = size(this.A,1);
      
      %% naive way
      %C = this.expand() * B;

      %% for-loop cell2mat
      %C = cell(this.d,1);
      %for ii = 0:(this.d-1)
      %  C{ii+1} = this.A * B(ii*n + (1:n),:);
      %end
      %C = cell2mat(C);

      %% for-loop arrays 
      %C = zeros(m*this.d,k);
      %for ii = 0:(this.d-1)
      %  C(ii*m + (1:m),:) = this.A * B(ii*n + (1:n),:);
      %end


      %% Transpose each "block" of B, so that if B=[B1;B2;...;Bd];
      %% now B = [B1 B2 ... Bd];
      %BB = reshape(permute(reshape(B',[k n this.d]),[1 3 2]),[this.d*k n]);
      %CC = BB * this.A';
      %% now we have made C = [A*B1 A*B2 ... A*Bd]
      %C = reshape(permute(reshape(CC,[k this.d n]),[1 3 2]),[k n*this.d])';


      %% Transpose each "block" of B, so that if B=[B1;B2;...;Bd];
      %% now B = [B1 B2 ... Bd];
      %BB = cell2mat(reshape(mat2cell(B,n*ones(1,this.d),k),1,this.d));
      %CC = this.A * BB;
      %% now we have made C = [A*B1 A*B2 ... A*Bd]
      %C = cell2mat(reshape(mat2cell(CC,m,k*ones(1,this.d)),this.d,1));

      %% cellfun
      %timesA = @(X) this.A * X;
      %BB = mat2cell(B,n*ones(1,this.d),k);
      %C = cell2mat(cellfun(timesA,BB,'UniformOutput',false));

      % Apply to columns of each block (disregarding order but then
      % reassembling)
      BB = reshape(B,n,k*this.d);
      CC = this.A*BB;
      C = reshape(CC,m*this.d,k);
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
