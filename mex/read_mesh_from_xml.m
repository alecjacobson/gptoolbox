% READ_MESH_FROM_XML Read a mesh from an xml file written using
% igl::serialize_xml V from the <vertices> tag and F from the <faces> tag. V
% may be written using CGAL::Epeck but will be converted summarily to
% double upon reading.
%
% [V,F] = read_mesh_from_xml(filename)
%
% Inputs:
%   filename  path to .xml file
% Outputs:
%   V  #V by 3 list of vertex positions
%   F  #F by 3 list of face indices into V
%
