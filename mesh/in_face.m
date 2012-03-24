function I = in_face(V,F,P)
  % IN_FACE test for each p in P whether it lies inside each f in F defined
  % over V
  % 
  % I = in_face(V,F,P)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim+1 list of element indices
  %   P  #P by dim list of query positions
  % Outputs:
  %   I  #P by #F matrix of bools
  %

  dim = size(V,2);
  assert(dim+1 == size(F,2));

  % number of elements 
  m = size(F,1);
  % number of query points 
  np = size(P,1);
  
  % triangle side lengths
  l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];

  B = zeros([np m dim+1]);
  for ii = 1:(dim+1)
    jj = mod(ii+1,dim+1)+1;
    kk = mod(ii,dim+1)+1;
    ljj = all_pairs_distances(P,V(F(:,jj),:));
    lkk = all_pairs_distances(P,V(F(:,kk),:));
    % semiperimeters
    s = bsxfun(@plus,l(:,ii)',ljj + lkk)*0.5;
    % Heron's formula for area
    B(:,:,ii) = 2*sqrt(s.*(bsxfun(@minus,s,l(:,ii)').*(s-ljj).*(s-lkk)));
  end
  % sum of barycentric coordinates
  sumA = sum(B,3);
  % area of element
  dblA = doublearea(V,F);
  %% check whether sum is more than true are
  %I = ~bsxfun(@gt,sumA,dblA');
  I = (bsxfun(@minus,sumA,dblA')) < sqrt(eps);
end

