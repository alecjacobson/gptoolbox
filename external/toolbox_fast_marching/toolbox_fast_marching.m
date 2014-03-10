% Toolbox Fast Marching - a toolbox for the computation of the Fast
% Marching algorithm both in 2D and 3D.
%
% The Fast Marching algorithm, introduced by Sethian (1996) is a numerical
% algorithm that is able to catch the viscosity solution of the Eikonal
% equation |grad(D)|=P. The level set {x \ F(x)=t} can be seen as
% a front advancing with speed P(x).
%
% The resulting function D is a distance function, and if the 
% speed P is constant, it can be seen as the distance function 
% to a set of starting points.
% 
% The Fast Marching is very similar to the Dijkstra algorithm
% that find shortes paths on graph. Using a gradient
% descent of the distance function D, one is able
% to extract a good approximation of the shortest path 
% (geodesic) in various settings (euclidean for P constant, 
% and a weighted riemanian manifold with P varying).
%
% The main reference about the Fast Marching algorithm is the book
%   Level Set Methods and Fast Marching Methods 
%   Evolving Interfaces in Computational Geometry, Fluid Mechanics, Computer Vision, and Materials Science
%   J.A. Sethian, Cambridge University Press, 1999
%   Cambridge Monograph on Applied and Computational Mathematics
%
% A good review of the Fast Marching in 3D together with some applications
% can be found in 
%   Fast extraction of minimal paths in 3D images and application to virtual endoscopy.  
%   T.Deschamps and L.D. Cohen. 
%   September 2000. To appear in  Medical Image Analysis.
%
% The function 'perform_fast_marching'
% computes the distance function from a set of starting points.
% To extract the geodesics between these starting points and 
% an ending point, you can use 'compute_geodesic'. 
%
% perform_fmstar_2d and perform_fmstar_3d implement the algorithm described in 
%	Heuristically Driven Front Propagation for Path Planning and Geodesic Extraction
%	Gabriel Peyré and Laurent Cohen, 2005 ,To appear.
% See script test_fmstar_2d and test_fmstar_3d for samples of use.
%
% The main computation are done in a mex file so it is very 
% fast (using a C++ heap structure). Precompiled version
% for windows are given.
%
%   Copyright (c) 2005 Gabriel Peyré