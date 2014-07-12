function [Nu,Nv] = neumannmatrix(mmesh,V, nlayers)
% NEUMANNMATRIX Computes matrices for computing neumann boundary conditions.
%
% [Nu,Nv] = neumannmatrix(mmesh,V)
%
% Inputs:
%   mmesh  hds representation of mesh
%   V  #V by dim list of mesh positions
%   nlayers  number of layers into mesh
% Outputs:
%   Nu  sparse matrices of size #V x #V,
%   Nv  sparse matrices of size #V x #V, The only nonzero columns/rows are for
%     boundary vertices.  The indices used in the matrix match the original
%     indices.  If the derivatives gradx_u = dx/du and gradx_v = dx/dv of the
%     unknown x are given at vertices of the boundary, we interpolate them
%     linearly and compute the vector withe entries <dx/dn phi_j>,  as \sum_i
%     <gradx dot t_ij, phi_j> = = \sum_i gradx_i dot t_ij <phi_i phi_j> = = (Nu
%     gradx_u + Nv gradx_v)_j, where phi_i are hats on the boundary and t_ij
%     are outward perpendiculars to the boundary edges
%
% See also: hds

    % indices of tails and tips of boundary edges
    % TODO should not be here
    layers = mmesh.layers; 
    interior = mmesh.interior; 
    layershe = mmesh.layershe;
    N0 = layers{nlayers};
    N0he = layershe{nlayers};
    tails = mmesh.tail(N0he)';
    tips  = mmesh.tip(N0he)';
    p1 = V(tails,:); p2 = V(tips,:);
    e = p2-p1;
    
    eperp = [e(:,2),-e(:,1)];
    % indices for nonzero entries 
    % per element 2x2 matrix, with indices tail, tip in each dimension
    i = [tails tips  tips tails];
    j = [tips  tails tips tails];
    % if (i,j) is a boundary edge with perp tij = [tij_x,tij_y], Nx(i,j) = tij_x/6
    % similarly Ny(i,j) = tij_y/6
    % and diag entries are  tij_x/3  + tik/3, where j and k are vertices
    % adjacent to i 
    val_u = [eperp(:,1)/6.0, eperp(:,1)/6.0, eperp(:,1)/3.0, eperp(:,1)/3.0];
    val_v = [eperp(:,2)/6.0, eperp(:,2)/6.0, eperp(:,2)/3.0, eperp(:,2)/3.0];    
    Nu = sparse(i,j,val_u, mmesh.nvert, mmesh.nvert);
    Nv = sparse(i,j,val_v, mmesh.nvert, mmesh.nvert);
end    
    
