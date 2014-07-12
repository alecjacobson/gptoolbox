function [polycoefs] = vtxpolyfit(V,mmesh,vtxset,vals,deg) 
% VTXPOLYFIT produces #V x m matrix of polynomial coefficients, m =
% (deg+1)*(deg+2)/2 the ordering of monomials is 1, u, v, u^2, uv v^2, u^3, u^2
% v ...  the fit is regularized, so that if there are not enough degrees of
% freedom a reasonable result should still be produced if 2 layers of vertices
% are fixed, there is enough dofs for deg 1 and typically not enough for deg 2
% if 3 layers are fixed then there is typically enough dofs for deg 2
%
% [polycoefs] = vtxpolyfit(V,mmesh,vtxset,vals,deg)
% 
% Inputs:
%   V  #V list of points
%   mmesh  data structure produced by hds
%   vtxset  a set of indices
%   vals  values at vertices
%   deg  polynomial total degree , primarily intended for deg = 1,2
% Outputs:
%   polycoefs  #V by m matrix of polynomial coefficients
% 

    num_monom = (deg+1)*(deg+2)/2;
    % powers of u and v for monomials u^l u^m, l+m <=deg, enumerated
    % sequentially for the least squares matrix
    powu = zeros(num_monom,1); powv = powu;
    last = 0;
    for i = 0:deg
        powu(last+1:last+i+1) = i:-1:0;
        powv(last+1:last+i+1) = 0:i;
        last = last+i+1;
    end;

    
    num_vtx = length(vtxset);
    vtxset_sorted = sort(vtxset);
    
    % **** build the lists of neighbors of each vertex in vtxset combined into one long vector

    % we only consider neighbors which are not in the interior
    noninterior = setdiff(1:mmesh.nvert,mmesh.interior);
    % adjacency matrix of vertices in vtxset restricted to outside of the
    % interior
    adj_vtxset = sparse([],[],[], mmesh.nvert,mmesh.nvert);
    adj_vtxset(vtxset_sorted,noninterior) = mmesh.adj(vtxset_sorted,noninterior);
    [vtx,nbrvtx,heind] = find(adj_vtxset);
    % extra vertices for fit stability: 
    % if he if the halfedge opposite a vertex 
    % we add the vertex of the triangle on the other side of the edge to 
    % the neighborhood, or an extra vertex along the boundary, if the 
    % edge is the boundary edge
    extravtx = mmesh.tip(mmesh.next(mmesh.opp(mmesh.prev(heind))));
    % can only use noninterior vertices
    validextras = ~ismember(extravtx,mmesh.interior);
    % create an adj matrix as if edges were connecting vertices to extras
    % so that it can be added to the adj_vtxset
    adj_extras = sparse( vtx(validextras), extravtx(validextras), 1, mmesh.nvert,mmesh.nvert); 
    %
    % add extra vertices
   adj_vtxset = adj_vtxset + adj_extras;
    % add the diagonal: each vertex is its own neighbor
    adj_vtxset = adj_vtxset + sparse(vtxset_sorted,vtxset_sorted,1,mmesh.nvert,mmesh.nvert);
    
    % nbrlists first col: vertex index; second col: neighbor vertex index
    [vtx,nbrvtx] = find(adj_vtxset);
    nbrlists = sortrows([vtx,nbrvtx],1);
    % nbrstarts: first indices of each 1-neighborhood list in nbrlists
    % vtxseqind:  the list of vertices in vtxset with each vertex index 
    % replaced by its sequential index in sorted nbrlists; 
    % (we need it to construct the monomials matrix)
    [dummy,nbrstarts, vtxseqind] = unique(nbrlists(:,1),'first');
    % ending indices of each neighbor list
    nbrends   = [nbrstarts(2:end); size(nbrlists,1)+1]-1;
    % number of neigbors for each vertex
    nbrsizes  = nbrends-nbrstarts+1;
    % computing it in a different way to verify correctness
    nbrsizes_alt = ( adj_vtxset ~= 0)*ones(size(mmesh.adj,2),1);
    nbrsizes_alt = nbrsizes_alt(nbrsizes_alt ~= 0);
    fprintf('verifying nbrhd sizes computed in a different way, diff should be zero: %g\n', ...
    norm(nbrsizes - nbrsizes_alt));
    if( any(nbrsizes < num_monom) ) 
      % TODO: add flaps, this will pretty much guarantee this does not
      % happen
       warning('vertex nbrhds have too few vertices');
       fprintf('%g ', find(nbrsizes < num_monom));
       fprintf('\n');
       fprintf('%g ', nbrsizes(nbrsizes < num_monom));
       fprintf('\n');
    end
        
    % ****** building a big matrix for all vertices at once
    
    % this is a block-diagonal matrix, with each diagonal block 
    % P^i of size #nbrs(i) x num_monom, i = 1.. num_vertices
    % row indices are just 
    nbrindices   = reshape(repmat(1:size(nbrlists,1)',num_monom,1),[],1);
    % build a matrix with first column being sequential indices of 
    % vertices in vtxset, each repeated as many times as there are vertices
    % in its neighborhood, followed by num_monom-1 columns of 1
    monomindices = [(vtxseqind-1)*num_monom+1 repmat(ones(size(vtxseqind,1),1),1,num_monom-1)];
    % compute cumm. sums  along rows: in this way, we get sequential indices 
    % in each row instead of 1; then flatten
    monomindices = reshape( cumsum(monomindices,2)',[],1);
    % finally, build the matrix monomial matrix
    P = sparse(nbrindices, monomindices, ...
        V(nbrlists(nbrindices,2),1).^powu(rem(monomindices-1,num_monom)+1).* ...
        V(nbrlists(nbrindices,2),2).^powv(rem(monomindices-1,num_monom)+1));
    % least squares system with basic regularization
    A = P'*P + 1e-15*speye(size(P,2));
    rhs = P'*vals(nbrlists(:,2));
    % least squares polynomial coefficients, arranged in one big vector, 
    % the order is the same as the sorted order of vertices in vtxset
    polycoefs = A \ rhs;
    % converting to a sparse matrix for ease of lookup
    polycoefs = sparse(reshape(repmat(vtxset_sorted,num_monom,1),[],1),repmat(1:num_monom,1,num_vtx),polycoefs,mmesh.nvert,num_monom);
end
    
