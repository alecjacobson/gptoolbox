function [OV,OF] = combine(IV,IF,EV,EF)
  % COMBINE Two meshes where the second mesh was constructed as the surface of
  % the union of the volume of the first mesh and a third mesh. The order of
  % vertices is such that the first chunk of vertices in the second mesh is a
  % subsequence of vertices in the first mesh (the remaining chunks are a
  % subsequence of vertices from the third mesh, then new vertices). The goal
  % is to create a new mesh where the vertices are: all the vertices from the
  % first mesh then other vertices and faces from the second mesh. This means
  % in the resulting mesh the unreferenced vertices correspond to vertices that
  % were in the first mesh but not in the second mesh.
  % 
  % Inputs:
  %   IV  vertex list of the first mesh
  %   IF  face list of the first mesh
  %   EV  vertex list of the second mesh
  %   EF  face list of the second mesh
  % Outputs:
  %   OV  vertex list of the output mesh
  %   OF  face list of the output mesh
  % 

  % create new list of vertices which is IV followed by vertices in EV not in
  % IV
  OV = [IV; setdiff(EV,IV,'rows')];
  % find location in new list for each vertex in EV
  [TF,LOC] = ismember(EV,OV,'rows');
  assert(min(LOC) > 0);
  OF = LOC(EF);
end
