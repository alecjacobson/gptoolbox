% ELTOPO Given a surface mesh (V0,F) in a self-intersecting state (may be many
% connected components) and desired new vertex positions (potentially causing
% self-intersections or "tunneling"), use EL TOPO's implementation of "ROBUST
% TOPOLOGICAL OPERATIONS FOR DYNAMIC EXPLICIT SURFACES" [Brochu & Bridson,
% 2009] to find new vertex positions U and possibly new combinatorics G
% resolving those collisions.
%
% [U,G,t] = eltopo(V0,F,V1)
%
% Inputs:
%   V0  #V by 3 list of input vertex positions at "time t=0"
%   F  #F by 3 list of triangle indices into V0
%   V1  #V by 3 list of desired output vertex positions at "time t=1"
% Outputs:
%   U  #U by 3 list of self-intersection free vertex positions at "time t=1"
%   G  #G by 3 list of triangle indices into U.
%   t  time step actually achieved (t=1 on success, t=0 on failure)
%
