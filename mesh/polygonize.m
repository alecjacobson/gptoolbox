function [CV,CE] = polygonize(V,F,fun)
  % Polygonize (contour) an implicit function in the spirit of "An Implicit
  % Surface Polygonizer" [Bloomenthal 1994]
  %
  % [CV,CE] = polygonize(V,F,fun)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions of original mesh
  %   F  #F by dim+1 list of element indices into V
  %   fun  function handle so that zero is the desired level set
  % Outputs:
  %   CV  #CV by dim list of contour mesh vertices 
  %   CE  #CE by dim list of facet indices into CV
  %
  max_iters = 10;
  D = fun(V);
  interval = @(DE) any(DE>0.0,2) & any(DE<=0.0,2);
  FF = F(interval(D(F)),:);
  % simplex size
  ss = size(F,2);
  switch ss
  case 3
    allE = [FF(:,[2 3]);FF(:,[3 1]);FF(:,[1 2])];
  case 4
    allE = ...
      [FF(:,1) FF(:,2); ...
       FF(:,1) FF(:,3); ...
       FF(:,1) FF(:,4); ...
       FF(:,2) FF(:,3); ...
       FF(:,2) FF(:,4); ...
       FF(:,3) FF(:,4) ...
       ];
  end
  [E,~,EMAP] = unique(sort(allE,2),'rows');
  crossing = interval(D(E));
  J = (1:size(E,1))';
  EE = E(crossing,:);
  J(crossing) = 1:size(EE,1);
  % Blasphemy
  switch ss
  case 3
    CE = sort(reshape(crossing(EMAP),[],3).*reshape(J(EMAP),[],3),2);
    CE = CE(:,2:end);
  case 4
    % CE(f,i) = 0 if ith edge of element f does not cross, otherwise
    % CE(f,i) is the index of the unique edge that does cross
    CE = reshape(crossing(EMAP),[],6).*reshape(J(EMAP),[],6);
    % If 3 edges cross then we can surface with a single triangle
    CT = sort(CE(sum(CE>0,2)==3,:),2);
    CT = CT(:,end-2:end);
    CQ = sort(CE(sum(CE>0,2)==4,:),2);
    CQ = CQ(:,end-3:end);
    CE = [CT;CQ(:,[1 4 3]);CQ(:,[1 2 4])];
  end
  assert(size(CE,2) == ss-1);
  % Upper and lower bound on barycenteric coordinate locating =0.5
  EEl= zeros(size(EE,1),1);
  EElV = V(EE(:,1),:);
  Dl = D(EE(:,1));
  EEu= ones(size(EE,1),1);
  Du = D(EE(:,2));
  EEuV = V(EE(:,2),:);
  for iter = 1:max_iters
    EEm = (EEl+EEu)/2;
    CV = 0.5*(EElV+EEuV);
    if iter < max_iters
      Dm = fun(CV);
      front = interval([Dl Dm]);
       EEu(front,:) =  EEm(front,:);
      EEuV(front,:) = CV(front,:);
       Du(front,:) =  Dm(front,:);
       EEl(~front,:) =  EEm(~front,:);
      EElV(~front,:) = CV(~front,:);
       Dl(~front,:) =  Dm(~front,:);
    end
  end
end
