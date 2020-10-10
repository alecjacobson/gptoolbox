% FORM_FACTOR Computer form factors for radiosity computation
%
% [F,V,BC,N] = form_factor(P,E)
%
% Inputs:
%   P  #P by 2 list of 2D vertices
%   E  #E by 2 list of edge indices into P
% Outputs:
%   F  #E by #E matrix of integrated form factor values
%   V  #E by #E matrix of visibility values
%   BC  #E by 2 list of edge barycenters
%   N  #E by 2 list of length-weighted normals
