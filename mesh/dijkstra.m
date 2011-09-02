function D = dijkstra( G , S )
% --------------------------------------------------------------------
%      Mark Steyvers, Stanford University, 12/19/00
% --------------------------------------------------------------------
%
% DIJKSTRA  Find shortest paths in graphs
% 	D = dijkstra( G , S ) use the full or sparse matrix G in which 
% 	an entry (i,j) represents the arc length between nodes i and j in a 
%	graph. In a full matrix, the value INF represents the absence of an arc; 
%	in a sparse matrix, no entry at (i,j) naturally represents no arc.
%		
%	S is the one-dimensional matrix of source nodes for which the shortest
%	to ALL other nodes in the graphs will be calculated. The output matrices
%	D and P contain the shortest distances and predecessor indices respectively.
%	An infinite distance is represented by INF. The predecessor indices contain
%	the node indices of the node along the shortest path before the destination
%	is reached. These indices are useful to construct the shortest path with the
%	function pred2path (by Michael G. Kay).
%
%	This function was implemented in C++. The source code can be compiled to a
%	Matlab compatible mex file by the command "mex -O dijkstra.cpp" at the Matlab
%	prompt. In this package, we provide a compiled .dll version that is 
%       compatible all Windows based machines.  If you are not working on a 
%       Windows platform, delete the .dll version provided and recompile from
%       the .cpp source file.  If you do not have the Matlab compiler or a Windows
%       platform, delete the .dll version and dijkstra will then call the
%	Matlab function dijk.m (by Michael G. Kay).  Note that this Matlab
%	code is several orders of magnitude slower than the C based mex file.

N = size( G , 1 );
D = dijk( G , S , 1:N );

