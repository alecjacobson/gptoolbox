# File naming convention #
The filename should be the same as the function name followed by `.m`. Use
lowercase with words separated by underscores: e.g. `adjacency_matrix.m`

## Exceptions ##
File I/O functions should be `readEXT` and `writeEXT` where `EXT is replaced
with the file extension. If the extension is ambiguous then it should be
`readEXT_program` and `writeEXT_program`, where `program` is replaced by the
short name of the program that generated that file.

# Function headers #
Each file contains a single exposed function and should contain a standard
header containing description of the function, the function prototype and
description of inputs/outputs. For example in `centroid.m`:

    function [C,vol] = centroid(V,F,varargin)
      % CENTROID Compute the centroid of a closed polyhedron bounded by (V,F)
      %
      % C = centroid(V,F)
      % [C,vol] = centroid(V,F,'ParameterName',ParameterValue, ...)
      %
      % Inputs:
      %   V  #V by 3 list of mesh vertex positions
      %   F  #F by 3 list of triangle mesh indices
      %   Optional:
      %     'Robust' followed by whether to use more robust but costlier method
      %       for nearly closed input. {false}
      % Outputs:
      %   C  3-vector of centroid location
      %   vol  total volume of polyhedron
      %
      ... % contents of function
    end

