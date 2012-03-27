function [ M] = massmatrix_intrinsic(l,F,nvert,masstype)
  % MASSMATRIX_INTRINSIC compute the mass matrix from edge lengths only
  %
  % [M] = massmatrix_intrinsic(l,F)
  %
  % Inputs:
  %  l: #F by 3, array of edge lengths of edges opposite each face in F
  %  F: #F by 3, list of indices of triangle corners
  %  nvert: number of vertices, only needed to set size
  %  masstype: full, barycentric, or voronoi
  % TODO: this is almost identical to massmatrix, 
  % only the area computation is different, need to refactor
  %
  % here's a handy line to view mass matrix entries on plot:
  % text(UV(:,1), UV(:,2),zeros(size(UV,1),1),num2str(M(M>0)))
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: massmatrix
  %

    % should change code below, so we don't need this transpose
    if(size(F,1) == 3)
      warning('F seems to be 3 by #F, it should be #F by 3');
    end
    F = F';

    % renaming indices of vertices of triangles for convenience
    l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
    % semiperimeters
    s = (l1 + l2 + l3)*0.5;
    % Heron's formula for area
    dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
    
    % renaming indices of vertices of triangles for convenience
    i1 = F(1,:); i2 = F(2,:); i3 = F(3,:); 
    
    if strcmp(masstype,'full')
        % arrays for matrix assembly using 'sparse'
        % indices and values of the element mass matrix entries in the order 
        % (1,2), (2,1),(2,3), (3,2), (3,1), (1,3) (1,1), (2,2), (3,3);
        i = [i1 i2 i2 i3 i3 i1  i1 i2 i3];
        j = [i2 i1 i3 i2 i1 i3  i1 i2 i3];
        offd_v = dblA/24.;
        diag_v = dblA/12.;
        v = [offd_v,offd_v, offd_v,offd_v, offd_v,offd_v, diag_v,diag_v,diag_v];  
    elseif strcmp(masstype,'barycentric')
        % only diagonal elements
        i = [i1 i2 i3];
        j = [i1 i2 i3];
        diag_v = dblA/6.;
        v = [diag_v,diag_v,diag_v];
    elseif strcmp(masstype,'voronoi')
      cosines = [ ...
      (l(:,3).^2+l(:,2).^2-l(:,1).^2)./(2*l(:,2).*l(:,3)), ...
        (l(:,1).^2+l(:,3).^2-l(:,2).^2)./(2*l(:,1).*l(:,3)), ...
        (l(:,1).^2+l(:,2).^2-l(:,3).^2)./(2*l(:,1).*l(:,2))];
      barycentric = cosines.*l;
      normalized_barycentric = barycentric./ ...
        [sum(barycentric')' sum(barycentric')' sum(barycentric')'];
      areas = 0.25*sqrt( ...
        (l(:,1) + l(:,2) - l(:,3)).* ...
        (l(:,1) - l(:,2) + l(:,3)).* ...
        (-l(:,1) + l(:,2) + l(:,3)).* ...
        (l(:,1) + l(:,2) + l(:,3)));
      partial_triangle_areas = normalized_barycentric.*[areas areas areas];
      quads = [ (partial_triangle_areas(:,2)+ partial_triangle_areas(:,3))*0.5 ...
        (partial_triangle_areas(:,1)+ partial_triangle_areas(:,3))*0.5 ...
        (partial_triangle_areas(:,1)+ partial_triangle_areas(:,2))*0.5];
    
      quads(cosines(:,1)<0,:) = [areas(cosines(:,1)<0,:)*0.5, ...
        areas(cosines(:,1)<0,:)*0.25, areas(cosines(:,1)<0,:)*0.25];
      quads(cosines(:,2)<0,:) = [areas(cosines(:,2)<0,:)*0.25, ...
        areas(cosines(:,2)<0,:)*0.5, areas(cosines(:,2)<0,:)*0.25];
      quads(cosines(:,3)<0,:) = [areas(cosines(:,3)<0,:)*0.25, ...
        areas(cosines(:,3)<0,:)*0.25, areas(cosines(:,3)<0,:)*0.5];

      i = [i1 i2 i3];
      j = [i1 i2 i3];
      v = reshape(quads,size(quads,1)*3,1);

    else 
        error('bad mass matrix type')
    end
    M = sparse(i,j,v,nvert, nvert);  
end
