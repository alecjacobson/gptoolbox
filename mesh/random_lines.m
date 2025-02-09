function [O,D,T,N] = random_lines(n,dim)
  % RANDOM_LINES Generate n lines uniformly randomly distributed in the
  % [-1,1]^dim hypercube
  %
  % [O,D] = random_lines(n,dim)
  %
  % Inputs:
  %   n  number of lines to generate
  %   dim  dimension of space in which to generate lines
  % Outputs:
  %   O  n by dim origin points on the surface of the hypercube
  %   D  n by dim direction vectors pointing into the hypercube
  %   T  n by 1 distance to other intersection with hypercube
  %

  D = normalizerow(randn(n,dim));
  
  % could special case this for dim=2
  switch dim
  case 2
    E = [D(:,2) -D(:,1)];
  case 3
    [~,~,V] = pagesvd(permute(D,[3 2 1]));
    E = permute(V(:,2:end,:),[3 1 2]);
  end

  %E = zeros(size(D,1),size(D,2),dim-1);
  %for i = 1:size(D,1)
  %  Di = D(i,:);
  %  Ei = null(Di);
  %  for j = 1:dim-1
  %    E(i,:,j) = Ei(:,j).';
  %  end
  %end
  %sum(D.*E,2)

  U = (rand(n,1,dim-1)*2-1)*sqrt(dim);

  O = sum(E.*U,3);

  O = O+D*sqrt(dim);

  T = nan;


  % determine if the ray O+D*t  intersects the [-1,1]^dim hypercube
  % using Slab method
  t_low = zeros(n,dim);
  t_high = zeros(n,dim);
  for i = 1:dim
    t_low(:,i) = (-1 - O(:,i))./D(:,i);
    t_high(:,i) = (1 - O(:,i))./D(:,i);
  end

  t_close = min(t_low,t_high);
  t_far = max(t_low,t_high);
  t_close = max(t_close,[],2);
  t_far = min(t_far,[],2);
  t = [t_close t_far];
  keep = t_close < t_far;
  if ~any(keep)
    O = [];
    D = [];
    T = [];
    N = [];
    return;
  end
  D = D(keep,:);
  O = O(keep,:);
  t_close = t_close(keep);
  t_far = t_far(keep);
  O = O + D.*t_close;
  T = t_far - t_close;
  N = size(O,1);

  if size(O,1) < n
    [Or,Dr,Tr,Nr] = random_lines(n-size(O,1),dim);
    O = [O;Or];
    D = [D;Dr];
    T = [T;Tr];
    N = [N;Nr];
  end

end
