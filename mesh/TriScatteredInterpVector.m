function [F] = TriScatteredInterpVector(V,D)
  % TRISCATTEREDINTERPVECTOR simple extension of TriScatteredInterp to handle
  % vector values at each data point
  %
  % [F] = TriScatteredInterpVector(V,D)
  %
  % Inputs:
  %   V  #V by 2 list of data positions
  %   D  #V by m list of data vectors
  % Outputs:
  %   A  m-size cell of TriScatteredInterp classes corrsepnding to each
  %     coordinate of data vectors
  %   ev  function handle to evaluate data at given points (X,Y)
  %

  m = size(D,2);

  % make room in cell array
  A = cell(1,m);
  % loop over data coordinates
  for ii = 1:m
    % build interpolant for this coordinate
    interp_func = TriScatteredInterp(V,D(:,ii));
    A{ii} = interp_func;
  end

  % function handle to eval interps
  F = @EvalTriScatteredInterpVector;

  function DXY = EvalTriScatteredInterpVector(X,Y)
    % for convenience build the interpolantion evaluation function handle
    f = @(I) I(X,Y);
    DXY = cellfun(f,A, 'UniformOutput', false);
    DXY = reshape(cell2mat(DXY),[size(X) m]);
  end
end
