function [VV,FF] = upsample(V,F,varargin)
  % UPSAMPLE Upsample a mesh by adding vertices on existing edges/faces
  %
  % [VV,FF] = upsample(V,F)
  % [VV,FF] = upsample(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #vertices by 3 list of vertex positions
  %   F  #faces by 3 list of face indices
  %   Optional:
  %     'KeepDuplicates' followed by either true or {false}}
  %
  % Outpus:
  %  VV new vertex positions
  %  FF new list of face indices
  %
  % This is Loop subdivision without moving the points
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  
%   % Add a new vertex at each face barycenter
%   % compute barycenters
%   C = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
%   % append barycenters to list of vertex positions
%   VV = [V;C];
%   % list of indices to barycenters
%   i = size(V,1) + (1:size(C,1))';
%   % New face indices, 3 new face for each original face
%   FF = [F(:,1) F(:,2) i; F(:,2) F(:,3) i; F(:,3) F(:,1) i];
%   

  keep_duplicates = false;
  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'KeepDuplicates'
      v = v+1;
      assert(v<=numel(varargin));
      keep_duplicates = varargin{v};
    otherwise
      error(['Unknown parameter: "' varargin{v} '"']);
    end
    v = v+1;
  end
  
  switch size(F,2)
  % http://mathoverflow.net/questions/28615/tetrahedron-splitting-subdivision
  case 3
    % Add a new vertex at the midpoint of each edge
    % compute midpoints (actually repeats, one midpoint per edge per face)
    m = [ ...
      (V(F(:,1),:) + V(F(:,2),:))/2; ...
      (V(F(:,2),:) + V(F(:,3),:))/2; ...
      (V(F(:,3),:) + V(F(:,1),:))/2];
    % indices of midpoints
    i3 = size(V,1) + (1:(size(m,1)/3))';
    i1 = (size(V,1)+(size(m,1)/3)) + (1:(size(m,1)/3))';
    i2 = (size(V,1)+(size(m,1)/3)+(size(m,1)/3)) + (1:(size(m,1)/3))';
    % new face indices, 4 new faces for each original face. As if we simply
    % ignored the duplicates in m and had appended m to V
    FF = [ F(:,1) i3 i2 ; F(:,2) i1 i3 ; F(:,3) i2 i1 ; i1 i2 i3];
    % find unique midpoints (potentially slow, though matlab is good at
    % these)
    if keep_duplicates
      m = m;
      J = (1:size(m,1))';
    else
      [m,~,J] = unique(m,'rows');
    end
    % append unique midpoints to vertex positions
    VV = [V ; m];
    % reindex map from duplicate midpoint indices to unique midpoint indices
    J = [(1:size(V,1))';J+size(V,1)];
    % reindex faces
    FF = J(FF);
  case 2
    m = [ (V(F(:,1),:) + V(F(:,2),:))/2 ];
    % indices of new midpoint vertices
    im = size(V,1) + (1:size(m,1))';
    % insert new face indices
    FF = [F(:,1) im;im F(:,2)];
    % append unique midpoints to vertex positions
    VV = [V;m];
    % No duplicates in 2D case
  end
  
end
