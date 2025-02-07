function [O,D,T] = random_lines(n,dim)
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
  [~,~,V] = pagesvd(permute(D,[3 2 1]));
  E = permute(V(:,2:end,:),[3 1 2]);

  %E = zeros(size(D,1),size(D,2),dim-1);
  %for i = 1:size(D,1)
  %  Di = D(i,:);
  %  Ei = null(Di);
  %  for j = 1:dim-1
  %    E(i,:,j) = Ei(:,j).';
  %  end
  %end
  %sum(D.*E,2)

  U = (rand(n,dim,dim-1)*2-1)*sqrt(dim);
  O = sum(E.*U,3);

  % determine if the ray O+D*t  intersects the [-1,1]^dim hypercube

  % find the t values for each dimension
  t = zeros(n,dim);
  for i = 1:dim
    t(:,i) = (1-O(:,i))./D(:,i);
    t(:,i+dim) = (-1-O(:,i))./D(:,i);
  end
  tpos = t;
  tpos(tpos<0) = inf;
  tpos = min(tpos,[],2);
  H_pos = O + D.*tpos;
  keep_pos = all(abs(H_pos) <= 1,2);

  tneg = t;
  tneg(tneg>0) = -inf;
  tneg = max(tneg,[],2);
  H_neg = O + D.*tneg;
  keep_neg = all(abs(H_neg) <= 1,2);

  keep = keep_pos | keep_neg;
  H = nan(n,dim);
  H(keep_pos,:) = H_pos(keep_pos,:);
  H(keep_neg,:) = H_neg(keep_neg,:);

  



  %[BV,BF] = bounding_box([-ones(n,dim);ones(n,dim)]);
  %clf;
  %hold on;
  %txt(O);
  %sct(O,'k','filled');
  %sct(O(keep,:),'b','filled');
  %sct(H(keep,:),'b','filled');
  %qvr(O, D,0,'k','LineWidth',1);
  %qvr(O,-D,0,'k','LineWidth',1);
  %qvr(O(keep,:), D(keep,:),0,'b','LineWidth',2);
  %qvr(O(keep,:),-D(keep,:),0,'b','LineWidth',2);
  %tsurf(BF,BV,'FaceColor','w','FaceAlpha',0.5);
  %hold off;
  %axis equal;

  O = H(keep,:);
  D = D(keep,:);
  N = (abs(O)==max(abs(O),[],2)).*O;
  D = -sign(sum(N.*D,2)).*D;

  %clf;
  %hold on;
  %sct(O,'k','filled');
  %qvr(O,D*2*sqrt(dim),0,'k','LineWidth',1);
  %tsurf(BF,BV,'FaceColor','w','FaceAlpha',0.5);
  %hold off;
  %axis equal;

  if size(O,1) < n
    [Or,Dr] = random_lines(n-size(O,1),dim);
    O = [O;Or];
    D = [D;Dr];
  end

  if nargout > 2
    N = (abs(O)==max(abs(O),[],2)).*O;
    T = zeros(n,dim);
    for i = 1:dim
      Pi = (1-O(:,i))./D(:,i);
      Pi(N(:,i)==1) = inf;
      Mi = (-1-O(:,i))./D(:,i);
      Mi(N(:,i)==-1) = inf;
      Mi(Mi<0) = inf;
      Pi(Pi<0) = inf;
      Ti = min(Pi,Mi);
      T(:,i) = Ti;
    end
    T = min(T,[],2);

    %clf;
    %hold on;
    %sct(O,'k','filled');
    %qvr(O,D*2*sqrt(dim),0,'k','LineWidth',1);
    %sct(O+D.*T,'r','filled');
    %tsurf(BF,BV,'FaceColor','w','FaceAlpha',0.5);
    %hold off;
    %axis equal;
  end

end
