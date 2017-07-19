function [bs,ts] = split_backfacing(ts)
  % SPLIT_BACKFACING Split the triangles of a trisurf mesh into front and back
  % facing based on the current view.
  %
  % Inputs:
  %   ts  list of trisurf handles
  % Outputs:
  %   bs  list of new back facing trisurf handles
  %   ts  list of front facing modified trisurf handles
  %
  bs = {};
  for ti = 1:numel(ts)
    t = ts{ti};
    BC = barycenter(t.Vertices,t.Faces);
    N = normals(t.Vertices,t.Faces);
    back = sum(N.*bsxfun(@minus,BC,campos),2)>0;
    b = copyobj(t,t.Parent);
    t.Faces = t.Faces(~back,:);
    b.Faces = b.Faces( back,:);
    bs{end+1} = b;
  end
end
