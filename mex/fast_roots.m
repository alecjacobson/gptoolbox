% FAST_ROOTS drop-in replacement for MATLAB's roots function, but vectorized and
% allows specifying range to consider.
%
% X = fast_roots(P, Xmin, Xmax)
%
% Inputs:
%   P  #P by deg+1 matrix of polynomial coefficients (in same order as MATLAB's
%   roots, polyval, etc. Higher degree terms first).
%   Xmin  #P|1 list of minimum x values to consider for each polynomial.
%   Xmax  #P|1 list of maximum x values to consider for each polynomial.
% Outputs:
%   X #P by deg matrix of roots for each polynomial, with NaNs for missing roots.
%
