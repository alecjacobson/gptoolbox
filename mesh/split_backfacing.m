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
      if size(t.FaceVertexCData,2) == 3
        Ctype = 'FaceVertexCData';
        C = t.FaceVertexCData;
      else
        Ctype = 'CData';
        C = t.CData;
      end
      set(t,'Faces',t.Faces(~back,:),Ctype,C(~back,:));
      set(b,'Faces',b.Faces( back,:),Ctype,C( back,:));
    else
      set(t,'Faces',t.Faces(~back,:));
      set(b,'Faces',b.Faces( back,:));
    end
    bs{end+1} = b;
  end
end
