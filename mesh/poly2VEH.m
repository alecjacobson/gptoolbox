function [V,E,H,h] = poly2VEH(poly)
  % POLY2VEH convert poly struct array to vertex, edge and hole lists
  %
  % [V,E,H,h] = poly2VEH(poly)
  % Input:
  %   poly struct array for each component of polygon boundary. Each struct
  %     contains a loop of vertices (x,y) and a flag for whether this is an
  %     inner hole
  % Output:
  %   V  #V by 2 list of polygon vertices
  %   E  #E by 2 list of polygon edge indices
  %   H  #H by 2 list of hole positions
  %
  % This has the problem that it is throwing away the information about which
  % sets of edges are inner boundaries
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: png2poly 
  %

  V = []; 
  E = [];
  H = [];
  H
  components = size(poly,2);
  % loop over components collecting vertices, edges, and holes
  for component_index = 1:components
    % component should constitute at least one triangle
    if(size(poly(component_index).x,1) >=3)
      component_vertices = ...
        [poly(component_index).x, poly(component_index).y];
      V = [V; component_vertices];
      component_E = ...
        [ 1:size(poly(component_index).x,1) ;...
        [size(poly(component_index).x,1), ...
        1:(size(poly(component_index).x,1)-1)]]';
      E = [E ; ...
        size(E,1) + component_E(:,1) , ...
        size(E,1) + component_E(:,2) ];
      if poly(component_index).hole == 1 
        % find point inside polygon to be hole marker
        poly(component_index)
        component_hole = ...
          point_inside_polygon(component_vertices);
        H = [H; component_hole];
      end
    end
  end

  if(size(H,1) >= size(poly,2))
    warning('Number of holes >= number of components. Something is wrong...');
  end
end
