% SEGMENT_GRAPH Given a sparse, square graph of edges weights segment the nodes
% of the graph into connected sub-components using the greedy merge-based
% method of "Graph Based Image Segmentation".
%
% C = segment_graph(A)
% C = segment_graph(A,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   A  #A by #A sparse, square matrix of edge weights
%   Optional:
%     'Threshold' followed by "C" threshold to use (paper writes that this
%       roughly corresponds to minimum size, though it's really just adding a
%       weight of size/C to components. In any case, increasing this will tend
%       to produce larger segments.
%     'MinSize' followed by the minimum size of an output component. This
%       constraint is enforced as a _post process_.
% Output:
%   C  #A by 1 list of component ids
%
% Example:
%   [V,F] = load_mesh('~/Dropbox/models/Cosmic blobs/Model9.off');
%   A = adjacency_dihedral_angle_matrix(V,F);
%   [AI,AJ,AV] = find(A);
%   A = sparse(AI,AJ,exp(abs(pi-abs(AV-pi))),size(A,1),size(A,2));
%   L = -(A - diag(sum(A,2)));
%   C = segment_graph(L,'Threshold',500,'MinSize',20);
%   tsurf(F,V,'CData',C);
%   colormap(cbrewer('Set1',(max(C))));
%   view(2);
%   axis equal;

