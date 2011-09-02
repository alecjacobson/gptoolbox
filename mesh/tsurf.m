function tsurf(F,V,vertex_indices,face_indices)
  % TSURF trisurf wrapper, easily plot triangle meshes with(out) face or vertex
  % indices
  %
  % tsurf(F,V,vertex_indices,face_indices)
  %
  % F: list of faces #F x 3
  % V: vertex positiosn #V x 3 or #V x 2
  % vertex_indices: show vertex indices on plot
  %                 0 -> off
  %                 1 -> text and grey background
  %                >1 -> text
  % face_indices: show face indices on plot
  %                 0 -> off
  %                 1 -> text and grey background
  %                >1 -> text
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: trisurf
  %


  if(~exist('vertex_indices'))
    % off by default
    vertex_indices = 0;
  end
  if(~exist('face_indices'))
    % off by default
    face_indices = 0;
  end

  % number of vertices
  n = size(V,1);

  % nuymber of dimensions
  dim = size(V,2);

  if(dim==2 || (dim ==3 && sum(abs(V(:,3))) == 0))
    V = [V(:,1) V(:,2) 0*V(:,1)];
    dim = 2;
  elseif(dim>3 || dim<2 ) 
    error('V must be #V x 3 or #V x 2');
    return;
  end

  trisurf(F,V(:,1),V(:,2),V(:,3));

  % if 2d then set to view (x,y) plane
  if( dim == 2)
    view(2);
  end

  if(face_indices==1)
    FC = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
    text(FC(:,1),FC(:,2),FC(:,3),num2str((1:size(F,1))'),'BackgroundColor',[.8 .8 .8]);
  elseif(face_indices)
    FC = (V(F(:,1),:)+V(F(:,2),:)+V(F(:,3),:))./3;
    text(FC(:,1),FC(:,2),FC(:,3),num2str((1:size(F,1))'));
  end

  if(vertex_indices==1)
    text(V(:,1),V(:,2),V(:,3),num2str((1:n)'),'BackgroundColor',[.8 .8 .8]);
  elseif(vertex_indices)
    text(V(:,1),V(:,2),V(:,3),num2str((1:n)'));
  end
  % uncomment these to switch to a better 3d surface viewing mode
  %axis equal; axis vis3d;
  axis(reshape([min(V(:,1:dim));max(V(:,1:dim))],1,dim*2));

end
