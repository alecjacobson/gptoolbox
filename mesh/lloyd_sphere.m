function [N,F] = lloyd_sphere(n)
  % LLOYD_SPHERE Construct a triangle mesh of the sphere with n reasonably well
  % distributed vertices.
  %
  % Inputs:
  %   n  number of points
  % Outputs:
  %   N  #N by 3 list of vertex positions
  %   F  #F by 3 list of face indices into N
  %

  %N = randsphere(n,'Method','trig');
  %N = normalizerow(N.^3);
  N = subdivided_sphere(ceil(log2(1/10*(10*(n)-20)^(1/2))));

  N = N(1:n,:);
  F = convhulln(N);
  if n <= 12
    % nothing more can be done
    return;
  end

  M = massmatrix(N,F,'voronoi');
  %t = tsurf(F,N);
  %caxis auto;
  %set(t,'CData',full(diag(M)));
  %colorbar;
  %axis equal;
  while true
    A = adjacency_matrix(F);
    A = A*M;
    %A = bsxfun(@rdivide,A,sum(A,2));
    A = spdiags (1./sum (A,2), 0, size(A,1), size(A,1)) * A ;
    N_prev = N;
    N = A*N;
    % subtract off center of mass  (needed for small n)
    N = bsxfun(@minus,N,diag(M)'*N./sum(diag(M)));
    N = normalizerow(N);
    F = convhulln(N);
    M = massmatrix(N,F,'voronoi');
    er = trace((N-N_prev)'*M*(N-N_prev));
    %set(t,'Vertices',N,'Faces',F, ...
    %  ... 'CData',full(diag(M)));
    %  'CData',full(sum(adjacency_matrix(F),2)),'FaceLighting','phong','FaceColor','interp');
    %axis equal;
    %drawnow;
    %title(sprintf('%g',er));
    if er < 1e-07
      break;
    end
  end
end
