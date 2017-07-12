% WIRE_MESH Construct a "wire" or "wireframe" or "strut" surface mesh, given a
% one-dimensional network of straight edges.
%
% [V,F] = wire_mesh(WV,WE)
% [V,F,J] = wire_mesh(WV,WE,'ParameterName',ParameterValue, ...)
% 
% Inputs:
%   WV  #WV by 3 list of vertex positions
%   WE  #WE by 2 list of edge indices into WV
%   Optional:
%   'Thickness' followed by diameter thickness of wire {0.1 average edge length}
%   'PolySize'  followed by number of sides on each wire (e.g., 4 would produce
%     wires by connecting rectangular prisms) {4}
% Outputs:
%   V  #V by 3 list of output vertices
%   F  #F by 3 list of output triangle indices into V
%   J  #F list of indices into [0,#WV+#WE) revealing "birth simplex" of
%     output faces J(j) < #WV means the face corresponds to the J(j)th
%     vertex in WV. J(j) >= #WV means the face corresponds to the
%     (J(j)-#WV)th edge in WE.
%
