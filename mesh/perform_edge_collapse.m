function [W,G,V2W] = perform_edge_collapse(V,F,C,t)
  % PERFORM_EDGE_COLLAPSE  Perform t edge collapses on mesh (V,F) based on
  % collapses C read from qslim log file
  %
  % [W,G] = perform_edge_collapse(V,F,C,t)
  %
  % Inputs:
  %   V  #V by 3 input mesh vertex positions
  %   F  #F by 3 input mesh triangle indices (1-indexed)
  %   C  #collapse list of edge collapse structs with fields
  %     .a index into V of vertex collapsed to
  %     .b index into V of vertex collapsed from
  %     .da  1 by 3 displacement of vertex a
  %     .rm list of indices in F of faces removed by this collapse
  %     .mod list of indices in F of faces modified by this collapse
  %   t  target number of faces {0 for all collapse}
  % Outputs:
  %   W  #W by 3 collapsed mesh vertex positions
  %   G  #G by 3 collapsed mesh triangle indices (1-indexed)
  %   V2W  #V list of indices into W
  %   
  % See also: readLOG, qslim
  %

  % collapsed versions
  G = F;
  W = V;
  rm_so_far = 0;
  if nargin < 4
    t = 0;
  end
  viz = false;
  if viz
    t = tsurf(G,W);
    axis equal;
    view(2);
  end
  V2W = (1:size(V,1))';
  % loop over collapses
  for ci = 1:numel(C)
    c = C(ci);
    % move vertex
    W(c.a,:) = W(c.a,:) + c.da;
    % remove faces in a cheapskate, in-place way
    G(c.rm,:) = 1;
    rm_so_far = rm_so_far + numel(c.rm);
    V2W(c.b) = c.a;
    % modify faces
    G(c.mod,:) = (G(c.mod,:) == c.b)*c.a + (G(c.mod,:) ~= c.b).*G(c.mod,:);
    if size(F,1)-rm_so_far <= t
      break;
    end
  
    if viz && (mod(ci,100)==0 || ci == numel(C))
      set(t,'Vertices',W,'Faces',G);
      title(sprintf('%d/%d',ci,numel(C)));
      drawnow;
    end
  end

  % peal off layers of indirection
  found = V2W==(1:size(V,1))';
  while any(~found)
    old_found = found;
    found(~old_found) = found(V2W(~old_found));
    V2W(~old_found) = V2W(V2W(~old_found));
  end
  
  G = G(~all(G==1,2),:);
  [W,IM] = remove_unreferenced(W,G);
  G = IM(G);
  V2W = IM(V2W);

end
