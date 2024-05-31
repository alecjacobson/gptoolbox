function [R,C,B] = circumradius(V,T)
  % CIRCUMRADIUS Return the circumradius of each triangle/tet element
  % 
  % R = circumradius(V,T)
  % [R,C,B] = circumradius(V,T)
  %
  % Input:
  %   V  #V by dim list of vertex positions
  %   T  #T by simplex-size list of tet indices
  % Output:
  %   R  #T by 1 list of simplex circumradii
  %   C  #T by dim list of simplex circumcenters
  %   B  #T by simplex-size list of barycentric coordinates so that: 
  %     C(t,:) = B(t,:) * V(T(t,:),:)
  %
  % Known issues:
  %   B output only supported for triangles
  %
  
  if size(T,2) == 3 && strcmp(class(V),'double')
    % http://www.mathworks.com/matlabcentral/fileexchange/17300-circumcircle-of-a-triangle/content/circumcircle.m

    %compute the length of sides (AB, BC and CA) of the triangle
    l = [ ...
      sqrt(sum((V(T(:,2),:)-V(T(:,3),:)).^2,2)) ...
      sqrt(sum((V(T(:,3),:)-V(T(:,1),:)).^2,2)) ...
      sqrt(sum((V(T(:,1),:)-V(T(:,2),:)).^2,2)) ...
      ];
    dblA = doublearea(V,T);
    %use formula: R=abc/(4*area) to compute the circum radius
    R = prod(l,2)./(2*dblA);
    %compute the barycentric coordinates of the circum center
    B= [ ...
      l(:,1).^2.*(-l(:,1).^2+l(:,2).^2+l(:,3).^2) ...
      l(:,2).^2.*( l(:,1).^2-l(:,2).^2+l(:,3).^2) ...
      l(:,3).^2.*( l(:,1).^2+l(:,2).^2-l(:,3).^2)];
    % normalize
    B = bsxfun(@rdivide,B,sum(B,2));
    %convert to the real coordinates
    C = zeros(size(T,1),size(V,2));
    for d = 1:size(V,2);
      C(:,d) = sum(B.*[V(T(:,1),d) V(T(:,2),d) V(T(:,3),d)],2);
    end
  else
    % https://math.stackexchange.com/a/4056112/35376
    % Augh, this was deleted _while_ I was looking at it.
    %
    % V ∈ ℝ^(d × (d+1)) corners
    % α ∈ ℝ^(d+1) barycentric coordinates of the 
    % C = Vα ∈ ℝ^d circumcenter
    % r ∈ ℝ circumradius
    % 
    % ‖vᵢ - Vα‖² = r² ∀ i = 1…d+1
    % ‖vᵢ‖² - 2 vᵢᵀVα + αVᵀVα = r² ∀ i = 1…d+1
    %
    % Let's call wᵢ = ‖vᵢ‖²
    %
    % And let λ = r² - αVᵀVα
    %
    % Then gather the equations into rows:
    %
    % 2 Vᵀ V α + λ = w
    %
    % We also have the constraint that ∑αᵢ = 1
    %
    % Together solve this system:
    %
    % [ 2VᵀV 𝟙 ] [ α ] = [ w ]
    % [ 𝟙ᵀ   0 ] [ λ ]   [ 1 ]
    %
    m = size(T,1);
    U = permute(reshape(V(T,:),size(T,1),size(T,2),size(V,2)),[3 2 1]);
    if m == 1
      UTU = U'*U;
    else
      UTU = pagemtimes(U,'transpose',U,'none');
    end
    A = 2*UTU;
    A(end+1,:,:) = 1;
    A(:,end+1,:) = 1;
    A(end,end,:) = 0;
    W = permute(sum(U.^2,1),[2 1 3]);
    b = W;
    b(end+1,:) = 1;
    if m == 1 && ~isnumeric(V)
      B = inv(A)*b;
    else
      if strcmp(class(V),'dlarray')
        B = my_pagemldivide(A,b);
      else
        B = pagemldivide(A,b);
      end
    end
    lambda = B(end,:,:);
    B = B(1:end-1,:,:);
    if m == 1
      C = permute(U*B,[3 1 2]);
    else
      C = permute(pagemtimes(U,B),[3 1 2]);
    end
    B = permute(B,[3 1 2]);
    R = sqrt(permute(lambda,[3 1 2]) + sum(C.*C,2));
  end

  function X = my_pagemldivide(A,b)
    n = size(A,1);
    solve_name = sprintf('solve%d',n);
    if ~exist(solve_name,'file')
      sA = sym('sA',[n,n,1],'real');
      sb = sym('sb',[n,1],'real');
      matlabFunction(simplify(inv(sA)*sb),'File',solve_name);
    end
    solve = str2func(solve_name);
    CA = mat2cell(A,ones(n,1),ones(1,n),size(A,3));
    Cb = mat2cell(b,ones(n,1),1,size(b,3));
    X = solve(CA{:},Cb{:});
  end
end
