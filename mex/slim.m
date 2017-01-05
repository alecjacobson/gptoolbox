% SLIM Comptue a parameterization using Scalable Locally Injective Maps
% [Rabinovich et al. 2016]
%
% U = slim(V,F,b,bc);
% U = slim(V,F,b,bc,'ParameterName',ParameterValue, ...)
%
% Inputs
%   V  #V by 3 list of 3D mesh coordinates
%   F  #F by 3 list of triangle indices into F
%   b  #b list of boundary indices into V
%   bc  #b by 2 list of boundary conditions
%   Optional:
%     'P'  followed by weight on boundary conditions: 0-->ignored,
%       large->enforced {1e5}
%     'Iters'  followed by number of iterations {100}
% Outputs:
%   U  #V by 2 list of output 2D mesh coordinates
%   U0  #V by 2 list of output 2D mesh coordinates of initial feasible guess
%
