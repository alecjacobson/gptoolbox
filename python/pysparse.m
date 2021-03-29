function outA = pysparse(inA)
  % Convert to and from python sparse matrix format.
  %
  % outA = pysparse(inA)
  % 
  % Inputs:
  %   inA  sparse matrix in either matlab or python format
  % Outputs:
  %   outA  if input is matlab, then python (CSC) sparse matrix; else if python
  %     then matlab sparse matrix.
  %
  if isa(inA,'double') && issparse(inA)
    [AI,AJ,AV] = find(inA);
    outA = py.scipy.sparse.csc_matrix( ...
      {AV,{uint64(AI-1) uint64(AJ-1)}}, ...
      {uint64(size(inA,1)),uint64(size(inA,2))});
    return;
  end
  if startsWith(class(inA),'py.scipy.sparse') && py.hasattr(inA,'tocoo')
    coo = inA.tocoo;
    outA = sparse(double(coo.row+1),double(coo.col+1),double(coo.data),double(coo.shape{1}),double(coo.shape{2}));
    return;
  end
  error('inA is unknown/non-sparse type (%s)',class(inA));
end
