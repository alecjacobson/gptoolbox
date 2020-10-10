function [k,H,K,M,T] = discrete_curvatures(V,F)
  % DISCRETE_CURVATURES Compute discrete mean, Gaussian and principal curvatures
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   k  #V by 2 list of integrated major and minor curvatures
  %   H  #V list of integrated mean-curvature values
  %   K  #V list of integrated Gaussian curvature
  %
  % See also:
  %   discrete_gaussian_curvature, discrete_mean_curvature
  %
  K = discrete_gaussian_curvature(V,F);
  H = discrete_mean_curvature(V,F);
  %a = -1;
  %b = 2*H;
  %c = -K;
  %k = (-b+[1 -1].*sqrt(b.^2 + 4*a*c))./(2.*a);
  %k2 = k;
  %k = k(:,1);
  %k = sort([k 2*H-k],2,'descend');
  M = massmatrix(V,F);
  %k = M*((M\H) + [1 -1].*sqrt((M\H).^2 - M\K));
  % Precision optimization:
  % h±m√((h/m)² - (k/m))
  % h±m√(h²/m² - km/m²)
  % h±m√(h²-km)/√m²
  % h±m√(h²-km)/m
  % h±√(h²-km)
  k = H + [1 -1].*sqrt(H.^2 - M*K);

  if nargout > 4
    % Copied from discrete_mean_curvature
    [A,C] = adjacency_dihedral_angle_matrix(V,F);
    [AI,AJ,AV] = find(A);
    [CI,CJ,CV] = find(C);
    assert(isequal(CI,AI));
    assert(isequal(CJ,AJ));
    % index into F(:)
    opp = sub2ind(size(l),CI,CV);
    inc = [ ...
      sub2ind(size(l),CI,mod(CV+1-1,3)+1) ...
      sub2ind(size(l),CI,mod(CV+2-1,3)+1)];
    lV = l(opp);

    T = zeros(size(V,1),3,2);
    for v = 1:size(V,1)
    end
  end

end
