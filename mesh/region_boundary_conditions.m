function [b,bc] = region_boundary_conditions(R,ii)
  % REGION_BOUNDARY_CONDITIONS compute boundary and boundary conditions for
  % solving for correspondence weights over a mesh (V,F) with control regions
  % defined on this mesh via ids in R, interior vertices (those not part of a
  % handle) will have id 0
  % 
  % [b,bc] = region_boundary_conditions(R,ii)
  % 
  % Inputs:
  %   R #V list of region ids, 0 meaning not part of any handle
  % Outputs:
  %   b  list of boundary vertices
  %   bc list of boundary conditions, size(boundary) by # handles matrix of
  %     boundary conditions where bc(:,i) are the boundary conditions for the 
  %     ith handle ( handle order is point handles then edges handles: P,E)
  %
  % Example:
  %   S = readSEL(sel_file_name);
  %   % convert survey style selection to regions
  %   R = (S==0)*1 + (S==2)*2;
  %   [b,bc] = region_boundary_conditions(R);
  %
  % See also: boundary_conditions, select_region, readSEL
  %

  assert(size(R,1) == 1 || size(R,2) == 1);
  % make list of ids into row
  R = R(:)';
  
  % number of mesh vertices
  n = numel(R);
  % indices to mesh vertices
  indices = 1:n;

  %b = [indices(R==0) indices(R==2)];
  b = indices(R~=0);
  %bc = zeros(numel(b),max(R(R~=0)));
  bc = sparse(1:numel(b),R(R~=0),1,numel(b),max(R(R~=0)));
end
