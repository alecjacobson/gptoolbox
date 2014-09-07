function [flag] = point_visible(V,F,o,heuristic_ratio)
  % POINT_VISIBLE  test whether vertices of mesh are visible to a set of 
  % points
  %
  % [flag] = point_visible(V,F,o)
  %
  % Input:
  %    V  #V by 3 list of vertex positions
  %    F  #F by 3 list of triangle indices
  %    o  row vector of position to test
  % Output:
  %    flag  #V by 1 list of bools (true) visible, (false) obstructed
  %
  
  dim = size(V,2);
  % Try to use accelerated bone_visible code
  if dim == 3 && 3==exist('bone_visible_embree','file')
    flag = bone_visible_embree(V,F,o,o);
  % otherwise try to use mex bone_visible code
  elseif dim ==3 && 3==exist('bone_visible','file')
    warning('Using non accelerated visibility test');
    flag = bone_visible(V,F,o,o);
  % otherwise use super slow matlab code
  else
    warning('Using non accelerated, pure matlab visibility test');
  
  
    assert(size(o,2) == dim);
    assert(size(F,2) == 3);
    % number of mesh vertices
    nv =  size(V,1);
  
    if dim == 2
      % get polygon edges of outline of 2D mesh
      O = outline(F);
    end
  
    if(exist('heuristic_ratio','var'))
      indices = randperm(nv);
      if(heuristic_ratio < 1)
        indices = indices(1:round(heuristic_ratio*nv));
      else
        indices = indices(1:round(heuristic_ratio));
      end
      flag = 0.5+eps + zeros(nv,1);
    else
      indices = 1:nv;
      flag = zeros(nv,1);
    end
  
    % loop over mesh vertices
    for jj_i = 1:numel(indices)
      progressbar(jj_i,numel(indices));
      jj = indices(jj_i);
      % query point, can o see q ?
      vjj = V(jj,:);
  
      dir = vjj-o;
  
  
      if dim == 3
        % extract triangles not containing jj
        Fmjj = F( all(F~=jj,2),:);
        %Fmjj = F;
        % shoot ray from o to q and see if any triangles not containing query get
        % hit
        %[hitc,tc] = mexray_mesh_intersect(o,dir,V,Fmjj);
        [hit,t] = ray_mesh_intersect(o,dir,V,Fmjj);
      elseif dim == 2
        Omjj = O( all(O~=jj,2),:);
        [hit,t] = ray_polygon_intersect(o,dir,V,Omjj);
      else
        error('Bad dimension');
      end
  
      %vjj-o
      %meshplot(V,F,'FC',repmat(hit,1,3))
      dir_mag =  sqrt(sum(dir.^2,2));
      flag(jj) = ~any(t(hit) < dir_mag);
    end
  
    if(exist('heuristic_ratio','var'))
      L = adjacency_matrix(F);
      L = L - diag(sum(L));
      i = 0;
      while(any(flag == (0.5+eps)))
        flag = flag+0.08*L*flag;
        flag(flag ~= (0.5+eps)) = round(flag(flag ~= (0.5+eps)));
        i = i +1;
      end
      i
      %L = cotmatrix(V,F);
      %I = speye(nv,nv);
      %L(indices,:) = I(indices,:);
      %flag = (L\flag)>0.5;
    end
  
  end
end
