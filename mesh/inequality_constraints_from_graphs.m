function [Aleq,ATleq] = inequality_constraints_from_graphs(F,A,AT)
  % INEQUALITY_CONSTRAINTS_FROM_GRAPHS  converts adjacency matrix graphs into
  % linear inequality constraints
  % 
  % [Aleq,ATleq] = inequality_constraints_from_graphs(F,A,AT)
  % 
  % Inputs:
  %  F  #triangles by 3 list of triangle indices
  %  A  #handles cell of #n by #n 
  %    monotonicity graph (sparse adjacency matrix) with A{k}(i,j) !=
  %    0 if we want that W(i,k) <= W(j,k), default is to use graph
  %    from harmonic coordinates
  %  AT #handles cell of #F by #F set of
  %    constraints on triangle (averages)
  % Outputs:
  %   Aleq  #nnz(A) by n  list of inequality constraints coming from A
  %   ATleq  #nnz(AT) by n  list of inequality constraints coming from AT
  %
  % See also: monotonic_biharmonic
  %
  assert(~isempty(A));
  n = size(A,1);

  [AI,AJ,AV] = find(A);
  nleq = size(AI,1);
  Aleq = ...
    sparse([1:nleq 1:nleq]',[AI;AJ],[ones(nleq,1);-ones(nleq,1)],nleq,n);
  ATleq = [];
  if nargin>=3 && ~isempty(AT)
    [ATI,ATJ,ATV] = find(AT);
    nleq = size(ATI,1);
    ATleq = ...
      sparse( ...
        repmat(1:nleq,1,6)', ...
        [F(ATI,1); F(ATI,2); F(ATI,3); F(ATJ,1); F(ATJ,2); F(ATJ,3);], ...
        [repmat( ones(nleq,1)/3,3,1); repmat(-ones(nleq,1)/3,3,1);], ...
        nleq, n);
  end
end
