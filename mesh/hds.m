% Halfedge datastructure from the list of triangles 
% mesh is assumed to be orientable, and triangle vertices enumerated 
% consistently counter-clockwise about a consistently chosen normal direction, 
% i.e, if two triangles share  the edge (i1,i2), the order of i1 and i2 in 
% the vertex list is (i1,i2) for one triangle and the opposite for the
% other.
%
% Halfege = edge + orientation; by convention, associated with the face
% for which traversal from tail to tip agrees with the order of vertices
% in the face list.
%
% All connectivity information is represented conceptually by 
% relationship between halfedges, and these are used as handles for 
% vertices and faces; i.e., to refer to a vertex v we use one of the 
% halfedges such that tip(he) = v, same for faces
%
% for each halfedge, we support the following const-time operations: 
% -- next halfedge around assoc. face
% -- prev halfedge around assoc. face
% -- opposite halfedge (i.e., the other halfedge for the same edge)
% -- tail and tip vertices 
% -- associated face
%
%prev(he1)\  he2face(he1)       /
%          \       he1         / next(he1)
%  tail(he1)*---------------->* tip(he1)
%           <----------------
%          /      opp(he1)   \
% 
%  border edges have special halfedges on one side, which do not have
%  the associated face.
%
%  this file implements the halfedge data structure for triangle meshes
%  only.
%

function mesh = hds(F,nlayers) 
% input: nf x 3 matrix of indices of vertices of mesh triangles 
% fields of the output struct (and corresponding MCGL fields)
% nhe     number of halfedges            ne
% tail   tail vertex of the halfedge     EO
% tip    tip  vertex of the halfedge     ED
% opp    opposite halfedge               ET
% he2face face corresp. to the halfedge  EF
% next   next halfedge around a face     EN
% prev   prev halfedge around a face     EP
% --                                     VE  can be done, for now do not need
% face2he                                FE  
% --                                     nb  (# boundary loops)
% --                                     BFS (size of each loop)
% --                                     IBF (bool flags for boundary faces)
% --                                     IBV (bool flags for boundary vert)
% bndhe array of indices of border       --
%        halfedges, in no particular order
% nfaces  number of faces                --
% nvert   number og vertices
% layers   cell array of nlayers (see getlayers)

nf = size(F,2); % number of faces
mesh.nfaces = nf;
mesh.F = F;
mesh.nhe = 3*nf;  % number of nonborder halfedges
mesh.face2he = [1:3:3*nf];  % first halfedge of each face
mesh.nvert = max(F(:));

% **** initialize all easy-to-do maps for nonborder halfedges ****
mesh.he2face = [repmat( 1:nf,1,3)]; % halfedge -> face table
mesh.prev = [(1:nf)+2*nf,(1:nf)     ,(1:nf)+nf];  
mesh.next = [(1:nf)+nf  ,(1:nf)+2*nf,(1:nf)];
mesh.tail = [ F(1,:) F(2,:) F(3,:)];
mesh.tip =  [ F(2,:) F(3,:) F(1,:)];

% **** now create the adjacency matrix and figure out the borders ***
he_indices = 1:mesh.nhe;
% he(i,j) = 1 if the (interior) halfedge exists
he_tag = sparse( mesh.tail, mesh.tip, ones(1,mesh.nhe));
% find halfedges on boundary edges 
[tip_b,tail_b,v_b] = find( he_tag - he_tag' == 1);
nbhe = size(tip_b,1);

% **** add border halfedges to easy-to-do maps
% add border halfedges (opposite the ones we've found) 
% to tail and tip functions
% their faces are set to zero
mesh.tail = [mesh.tail tail_b'];
mesh.tip  = [mesh.tip tip_b'];
mesh.he2face = [mesh.he2face, zeros(1,nbhe)];
bhe_indices = (mesh.nhe+1):mesh.nhe+nbhe;
mesh.nhe = mesh.nhe+nbhe;
mesh.bndhe = bhe_indices;


% prev and next for border halfedges (which do not have faces) 
% are a bit trickier: need to find the next/prev border halfedge
% along the boundary
% map tip -> border he index
tip2bhe =  sortrows([ tip_b,bhe_indices'],1);
% map tail -> border he index
tail2bhe = sortrows([ tail_b,bhe_indices'],1);  
% glue together, for i-th row, tip(blink(i,1)) = tail(blink(i,2))
blink = [tip2bhe(:,2), tail2bhe(:,2)];
next_b = sortrows( blink , 1); next_b = next_b(:,2);
prev_b = sortrows( blink , 2); prev_b = prev_b(:,1);
mesh.next = [mesh.next next_b'];
mesh.prev = [mesh.prev prev_b'];

% **** now build the opposite halfedge map
% he(i,j) = index of halfedge (i,j)
he = sparse( mesh.tail, mesh.tip, [he_indices bhe_indices]);
% after transpose adjmat(i,j) = index of halfedge (j,i) 
mesh.adj = he';
% find creates triple representations of matrices
[he_i  he_j  he_v ] = find(he);
[adj_i adj_j adj_v] = find(mesh.adj);
% after sorting, the first two columns of he and adjt are identical 
% and correspond to tails/tips of edges and the third column is the 
% index of the halfedge/opposite halfedge
het  = sortrows( [he_i he_j he_v], [1,2]);
adjt = sortrows( [adj_i adj_j adj_v], [1,2]);
% stitch together to create halfedge -> opposite map
% sort on the first row  
mesh.opp = sortrows( [het(:,3), adjt(:,3)],1);
mesh.opp = mesh.opp(:,2);
% 
% collect vertices and halfedges for layers (1 = boundary, 2 = 
% 1 edge away from the boundary etc, see getlayers
[mesh.layers,mesh.interior,mesh.layershe] = getlayers(mesh,nlayers);
end