function [F,Q,was_flipped] = shor_van_wyk(U)
  % By construction, determine whether a closed polyline path can be
  % triangulated using the Shor and van Wyk algorithm for self-overlapping
  % triangulations. "Detecting and decomposing self-overlapping curves" [Shor
  % and Van Wyk, 1991].
  % 
  % Inputs:
  %   U  #U by 2 list of points on a (possibly self-overlapping) closed polygon
  % Outputs:
  %   F  #U-2|0 by 3 list of faces (triangles) in the triangulation of the
  %     polygon), 0 if the polygon cannot be triangulated
  %   Q  #U by U shor and van wyk dynamic programming table.
  %   was_flipped  boolean, true if the polygon is negatively oriented. This is
  %     numerically unreliable if the polygon cannot be triangulated. (e.g., if
  %     the turning number is 0)
  % 
  % Example:
  %   [F,~,was_flipped] = shor_van_wyk(U);
  %   if was_flipped
  %      F = fliplr(F);
  %   end
  %   [dV,dF] = delaunayize(U,F);
  %   if was_flipped
  %      dF = fliplr(dF);
  %      F = fliplr(F);
  %   end
  %   

  function res = valid_dangle(j,i,k,d)
    res = ~(orient2d(U(j,:),U(i,:),U(d,:))>0 && orient2d(U(k,:),U(j,:),U(d,:))>0);
  end

  function helper(i,j,tab)
    %fprintf('%shelper(%d,%d)\n',tab,i,j);
    if ~isnan(Q(i,j))
      %fprintf('%s  we already know that Q(%d,%d) = %d\n',tab,i,j,Q(i,j));
      return;
    end
    Q(i,j) = 0;
    k = i;
    while true
      k = inc(k,1);
      if k == j
        break;
      end
      %fprintf('  Q(%d,%d) @ %d = %d\n',i,j,k,Q(i,j));

      if orient2d(U(j,:),U(i,:),U(k,:)) <= 0
        continue;
      end
      if ~valid_dangle(j,i,k,inc(j,-1)) ...
        || ~valid_dangle(i,k,j,inc(i,1)) ...
        || ~valid_dangle(k,j,i,inc(k,-1)) ...
        || ~valid_dangle(k,j,i,inc(k,+1))
        continue;
      end
      % other line segment tests...
      helper(i,k,[tab '  ']);
      helper(k,j,[tab '  ']);

      %fprintf('%s Q(%d,%d) = %d\n',tab,i,k,Q(i,k));
      %fprintf('%s Q(%d,%d) = %d\n',tab,k,j,Q(k,j));
      if Q(i,k) && Q(k,j)
        Q(i,j) = k;
        %fprintf('%s  Q(%d,%d) = %d !\n',tab,i,j,Q(i,j));
        %plt(U([1:end 1],:),'-or');
        %hold on;
        %sct(U(j,:),'filled');
        %sct(U(i,:),'filled');
        %sct(U(k,:),'filled');
        %tsurf([1 2 3],U([j i k],:));
        %txt(U([j i k],:),{'j','i','k'},'color','k','BackgroundColor','w','FontSize',14);
        %hold off;
        %axis equal;
        %pause
        break;
      end
    end
  end

  function [F,count] = build_F(i,j,inc,Q,F,count)
    k = Q(i,j);
    assert(k ~= 0, 'Q(i,j) must be non-zero');
    if k == -1
      return;
    end
    count = count+1;
    F(count,:) = [j i k];
    assert(Q(i,k) ~= 0 && Q(k,j) ~= 0, 'Q(i,k) and Q(k,j) must be non-zero');
    [F,count] = build_F(i,k,inc,Q,F,count);
    [F,count] = build_F(k,j,inc,Q,F,count);
  end


  %tic;
  n = size(U,1);
  I = 1:n;
  inc = @(I,a) mod(I+a-1,n)+1;
  was_flipped = false;
  [~,alpha] = curvature(U);
  if sum(alpha)<0
    U = flipud(U);
    was_flipped = true;
  end

  Q = nan(n,n);
  % Q(i,j) = 1 iff it is possible to triangule vᵢ,…,vⱼ
  % Q(i,i+1) = 1 (as a vacuous base case)
  % Q(i,i+2) = 1 iff vᵢ,vᵢ₊₁, vᵢ₊₂ form a convex angle at vᵢ₊₁
  % Since Q(i-1,1) = 1, then
  % Q(i,i-1) would mean that you can triangulate the other part of the polygon
  Q(sub2ind(size(Q),I,inc(I,1))) = -1;
  %[~,alpha] = curvature(U);
  alpha = orient2d(U(inc(I,-1),:),U(I,:),U(inc(I,1),:));
  positive = alpha(inc(I,1))>0;
  J = sub2ind(size(Q),I,inc(I,2));
  Q(J(~positive)) = 0;
  Q(J(positive)) = inc(I(positive),1);
  for i = 1:n
    %fprintf('i: %d\n',i);
    j = inc(i,-1);
    helper(i,j,'  ');
    if Q(i,j)
      break;
    end
  end
  %fprintf('%20s: %g secs\n','shor_van_wyk',toc);
  if Q(i,j) == 0
    %fprintf('sadly, we cannot triangulate the polygon\n');
    F = [];
  else
    %tic;
    F = zeros(n-2,3);
    count = 0;
    F = build_F(i,j,inc,Q,F,count);
    assert(size(F,1) == n-2, 'F must have n-2 rows');
    %fprintf('%20s: %g secs\n','build_F_transpose',toc);
    if was_flipped
      F = n-fliplr(F)+1;
    end
  end

end
