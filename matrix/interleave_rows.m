function C = interleave_rows(A,B)
  % Takes two nxm matrices and creates a 2nxm matrix such that even rows are
  % from the first matrix and odd rows are from the second.
  %
  % C = interleave_rows(A,B)
  %
  % Input:
  %   A  First nxm matrix, top row of output will be same as top of this matrix
  %      along with every other row
  %   B  Second nxm matrix, bottom row of output will be same as bottom of this
  %      matrix along with every other row
  %
  % Output:
  %   C  2nxm matrix, alternating rows
  %
  %
  if(size(A)~=size(B))
    error('Input matrices must be same shape');
  end
  C = reshape([A(:) B(:)]',size(A,1)*2, size(A,2));
end
