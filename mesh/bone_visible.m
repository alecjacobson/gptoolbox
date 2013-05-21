function [flag] = bone_visible(V,F,s,d);
  % BONE_VISIBLE  test whether vertices of mesh are visible to a set of "bones"
  % aka line segments. Here visibility is defined as [Baran and Popovic 2007]
  % defined it. A point q on the mesh is seen by bone [s,d] if the line segment
  % between q and the projection of q on [s,d] is entirely inside the mesh.
  %
  % [flag] = bone_visible(V,F,s,d);
  %
  % Input:
  %    V  #V by 3 list of vertex positions
  %    F  #F by 3 list of triangle indices
  %    s  row vector of position of start end point of bone
  %    d  row vector of position of dest end point of bone
  % Output:
  %    flag  #V by 1 list of bools (true) visible, (false) obstructed
  %
  

  dim = size(V,2);
  assert(size(s,2) == dim);
  assert(size(d,2) == dim);
  assert(size(F,2) == 3);
  if dim == 3 && 3==exist('bone_visible_embree','file')
    bone_visible_embree(V,F,s,d);
    return;
  elseif dim == 3 && 3==exist('bone_visible_embree','file')
    warning('Using non accelerated visibility test');
    bone_visible_mex(V,F,s,d);
    return;
  else
    warning('Using non accelerated, super slow pur matlab visibility test');

    % number of mesh vertices
    nv =  size(V,1);

    if dim == 2
      % get polygon edges of outline of 2D mesh
      O = outline(F);
    end

    flag = zeros(nv,1);
    % loop over mesh vertices
    for jj = 1:nv
      progressbar(jj,nv);
      % query point, can o see q ?
      vjj = V(jj,:);
      % project to line segment to find "closest origin" of bone
      t = project_to_lines(vjj,s,d);
      % snap to line segment
      t = min(max(t,0),1);
      % lerp to get position on line segment
      o = (1-t).*s+t.*d;
      dir = vjj-o;
      if dim == 3
        % extract triangles not containing jj
        Fmjj = F( all(F~=jj,2),:);
        %Fmjj = F;
        % shoot ray from o to q and see if any triangles not containing query get
        % hit
        [hit,t] = ray_mesh_intersect(o,dir,V,Fmjj);
      elseif dim == 2
        Omjj = O( all(O~=jj,2),:);
        [hit,t] = ray_polygon_intersect(o,dir,V,Omjj);
      else
        error('Bad dimension');
      end
      %vjj-o
      %hitF = 0;
      %hitF(all(F~=jj,2),:) = hit;
      %meshplot(V,F,'FC',repmat(hitF,1,3))
      %meshplot(V,F,'FC',repmat(hit,1,3))
      dir_mag =  sqrt(sum(dir.^2,2));
      flag(jj) = ~any(t(hit) < dir_mag);
    end
  end

end
