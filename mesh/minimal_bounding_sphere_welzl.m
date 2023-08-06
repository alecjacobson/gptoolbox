function [c,r] = minimal_bounding_sphere_welzl(P)
  % [c,r] = minimal_bounding_sphere_welzl(P)
  %
  % Inputs:
  %   P  #P by dim list of points
  % Outputs:
  %   c  1 by dim center of minimal bounding sphere
  %   r radius of minimal bounding sphere
  %

  % https://en.wikipedia.org/wiki/Smallest-circle_problem#Welzl's_algorithm

  function [c,r] = trivial(P,RI)
    assert(numel(RI) <= size(P,2)+1);
    switch numel(RI)
    case 2
      c = 0.5*(P(RI(1),:)+P(RI(2),:));
      r = 0.5*norm(P(RI(1),:)-P(RI(2),:));
    case 1
      c = P(RI,:);
      r = 0;
    case 0
      c = zeros(1,size(P,2));
      r = 0;
    otherwise
      % circumcircle
      [r,c] = circumradius(P,RI);
    end
  end

  function [c,r] = minimal_bounding_sphere_welzl_helper(P,k,RI)
    % Base case
    if k == 0 || numel(RI) == size(P,2)+1
      [c,r] = trivial(P,RI);
      return;
    end
    [c,r] = minimal_bounding_sphere_welzl_helper(P,k-1,RI);

    if sum((P(k,:)-c).^2,2) <= r*r
      return;
    end
    [c,r] = minimal_bounding_sphere_welzl_helper(P,k-1,[RI k]);
  end

  [c,r] = minimal_bounding_sphere_welzl_helper(P(randperm(size(P,1)),:),size(P,1),[]);
end
