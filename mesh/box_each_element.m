function [B1,B2] = box_each_element(V,F)
  % BOX_EACH_ELEMENT
  %
  % [B1,B2] = box_each_element(V,F)
  %
  % Inputs:
  %   V by dim list of mesh vertex positions
  %   F by ss list of simplices indexing rows of V
  % Outputs:
  %   B1 by dim list of min corners of boxes
  %   B2 by dim list of max corners of boxes
  %
  % corners of Axis-aligned boxes containing each element
  [B1,B2] = bounds(reshape(V(F,:)',size(V,2),[],size(F,2)),3);
  B1 = B1';
  B2 = B2';
  
end
