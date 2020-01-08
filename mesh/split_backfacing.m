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
    face_color = (size(t.FaceVertexCData,1) == size(t.Faces,1) || size(t.CData,1) == size(t.Faces,1)) && strcmp(t.FaceColor,'flat');
    if face_color
      set(t,'Faces',t.Faces(~back,:),'FaceVertexCData',t.FaceVertexCData(~back,:),'CData',t.CData(~back,:));
      set(b,'Faces',b.Faces( back,:),'FaceVertexCData',b.FaceVertexCData( back,:),'CData',b.CData( back,:));
    else
      set(t,'Faces',t.Faces(~back,:));
      set(b,'Faces',b.Faces( back,:));
    end
    bs{end+1} = b;
  end
end
