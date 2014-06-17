function Q = quadrature_points(V,F,k)
  % QUADRATURE_POINTS Sample quadrature points on a triangle or tetrahedron.
  %
  % Q = quadrature_points(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by simplex-size list of indices into V
  %   k  number of quadrature points {1,3,6}
  % Outputs:
  %   Q  #F by dim by k list of quadrature points
  %   w  #k list of weights
  %
  % Example:
  %   Q = quadrature_points(V,F);
  %   QQ = reshape(permute(Q,[1 3 2]),[],3);
  %   QQQ = permute(reshape(QQ,permute(size(Q),[1 3 2])),[1 3 2]);
  %   max(abs(Q-QQQ))
  % 
  % See also: barycenter
  %

  % Compute element quadrature points
  % "A SET OF SYMMETRIC QUADRATURE RULES ON TRIANGLES AND TETRAHEDRA", Table 2.2
  dim = size(V,2);
  ss = size(F,2);
  switch ss
  case 3
    switch k
    case 1
      lambda = [1/3 1/3 1/3];
    case 3
      %a = 0.5;
      a = 0.15;
      lambda = unique(perms([a a 1-2*a]),'rows');
    case 6
      a = 0.05;
      b = 0.1;
      lambda = unique(perms([a b 1-a-b]),'rows');
    otherwise
    error(sprintf('Unsupported number of quadratures %d',k));
    end
    % uniform weights
    w = ones(size(lambda,1),1)/size(lambda,1);
    QQ = cat(4,V(F(:,1),:),V(F(:,2),:),V(F(:,3),:));
  otherwise
    error(sprintf('Unsupported simplex-size %d',ss));
  end
  lambda = permute(lambda,[3 4 1 2]);
  Q = sum(bsxfun(@times,QQ,lambda),4);
end
