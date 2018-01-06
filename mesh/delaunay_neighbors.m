function N = delaunay_neighbors(q,P)
  % DELAUNAY_NEIGHBORS Find the Delaunay neighbors of q in a list of points P.
  % 
  % N = delaunay_neighbors(q,P)
  % 
  % Experimentally for uniformly randomly distributed points the number of
  % iterations is O(|P|) and specifically seems to do 1.355*|P| "iterations".
  % There's a sort on the angles so at best this implementation is O(|P|log|P|),
  % but it's very simple... so maybe the constants and overheads are much
  % smaller than doing a full Delaunay triangulation.
  %
  % Inputs:
  %   q  2D position of query point
  %   P  #P by 2 list of nearby points
  % Outputs:
  %   N  #N list of indices into P of neighboring points
  %

  % Example:
  %   q = [0 0];
  %   P = rand(60,2)*2-1;
  %   N = delaunay_neighbors(q,P);
  %   Nm = find(adjacency_matrix(delaunay([P;q]))*sparse(size(P,1)+1,1,1,size(P,1)+1,1));
  %   assert(isempty(setxor(N,Nm)));

  
  % 2----1
  % |θ₂ / \
  % |  /   \
  % | /     \
  % |/θq   θ₀\
  % q---------0
  %% Angle at p0-q-p1
  %thq = th(p1)-th(p0);
  %% Law of sines: θq/|01| = θ₀/|q1|
  %% θ₀ = |q1|*θq/|01|
  %th0 = l(p1)*thq/sqrt(sum((P(p0,:)-P(p1,:)).^2,2));
  % Vector from p0 to p1

  % Inputs:
  %   P  #P by 2 list of nearby points
  %   l  #P by 2 list of distances to origin: l = normrow(P)
  %   p0  first 2d point
  %   p1  second 2d point
  %   p2  third 2d point
  % Outputs:
  %   flag  whether edge from origin to p1 forms a Delaunay edge between points
  %     p1 and p2.
  % 
  isdel = @(P,l,p0,p1,p2) ...
    acos((-P(p0,:)/l(p0))*(((P(p1,:)-P(p0,:))/sqrt(sum((P(p1,:)-P(p0,:)).^2,2)))')) + ...
    acos((-P(p2,:)/l(p2))*(((P(p1,:)-P(p2,:))/sqrt(sum((P(p1,:)-P(p2,:)).^2,2)))')) <= pi;

  %% Or perhaps: http://stackoverflow.com/a/8523979
  %function f = isdel(P,l,p0,p1,p2)
  %  f = acos((-P(p0,:)/l(p0))*(((P(p1,:)-P(p0,:))/sqrt(sum((P(p1,:)-P(p0,:)).^2,2)))')) + ...
  %  acos((-P(p2,:)/l(p2))*(((P(p1,:)-P(p2,:))/sqrt(sum((P(p1,:)-P(p2,:)).^2,2)))')) <= pi;
  %  % Shewchuk's incircle test via determinant
  %  h(1,:) = [P(p0,:) 1];
  %  h(2,:) = [P(p1,:) 1];
  %  h(3,:) = [P(p2,:) 1];
  %  d = det([ ...
  %    q(1) q(2) q(1)^2+q(2)^2 1 ; ...
  %    h(1,1) h(1,2) h(1,1)^2+h(1,2)^2 h(1,3) ; ...
  %    h(2,1) h(2,2) h(2,1)^2+h(2,2)^2 h(2,3) ; ...
  %    h(3,1) h(3,2) h(3,1)^2+h(3,2)^2 h(3,3) ; ...
  %    ]);
  %  if f ~= (d<0)
  %    fprintf('%d %d\n',f,d<0);
  %    pause
  %  end
  %  %if any(l([p0 p1 p2]) > 1e3)
  %  %%if all((l([p0 p1 p2])>1e3)==[0;1;0]) && f
  %  %  fprintf('%d %d %d --> %d\n',l([p0 p1 p2]) > 1e5,f);

  %  %  clf;
  %  %  hold on;
  %  %  tsurf([1 2 3;1 3 4],[0 0;P([p0 p1 p2],:)],'FaceColor','none','LineWidth',2,'EdgeColor','b');
  %  %  tsurf([1 2 4;1 4 2],[0 0;P([p0 p1 p2],:)],'FaceColor','none','LineWidth',1,'EdgeColor','r');
  %  %  scatter([0;P([p0 p1 p2],1)],[0;P([p0 p1 p2],2)],'.','SizeData',500);
  %  %  text([0;P([p0 p1 p2],1)],[0;P([p0 p1 p2],2)],['q';num2str((0:2)')],'FontSize',20);
  %  %  hold off;
  %  %  axis equal;
  %  %  axis([-3 3 -3 3]);
  %  %  drawnow;

  %  %  pause
  %  %end
  %end

  % Subtract q from P and q
  P = bsxfun(@minus,P,q);
  q = [0 0];
  % Compute distances to origin
  l = sqrt(sum(P.^2,2));

  % Handle case where q is on the convex hull. 
  % TODO: Shouldn't _always_ have to added bounding points.
  max_l = max(l);
  % This is HUGE number is a hack. Should handle boundary "protectors" explicitly then.
  % Something like angle between:
  %   finite-∞-finite --> 0 
  %   ∞-finite-∞ --> pi
  %   finite-finite-∞ --> ?
  % Maybe this is a case for homogenous coordinates?
  s = 1e10*max_l;
  extra = [size(P,1)+(1:4)];
  P = [P;s*[1 1;1 -1;-1 -1;-1 1]];
  l = [l;repmat(s*sqrt(2),4,1)];
  assert(size(P,1)==size(l,1));

  % Number of input points
  n = size(P,1);
  % Compute angle with x-axis
  th = atan2(P(:,2),P(:,1));
  [th,I] = sort(th);
  [~,closest] = min(l);

  % Reorder so that closest comes first
  I = I([find(I==closest):end 1:find(I==closest)-1]);
  % Assumption: closest point must be Delaunay neighbor
  N = [I(1)];
  % Consider next two points
  p1i = 2;
  p2i = 3;
  p1 = I(p1i);
  p2 = I(p2i);
  p0 = N(end);
  k = 0;
  while true
    k = k+1;
    % is the edge q-p1 between p0 and p2 is delaunay
    if isdel(P,l,p0,p1,p2)
      % push p1 onto "keepers stack" N
      N = [N;p1];
    else
      % edge q-p1 between p0 and p2 is **not** delaunay
      % first try to work backward until finding delaunay
      while numel(N)>1 && ~isdel(P,l,N(end-1),N(end),p2)
        % MOVE BACKWARD
        k = k+1;
        p0 = N(end-1);
        p1 = N(end);
        N = N(1:end-1);
      end
    end
    if p1i == n
      break;
    end
    % MOVE FORWARD
    p1i = mod(p1i,n)+1;
    p2i = mod(p2i,n)+1;
    p1 = I(p1i);
    p2 = I(p2i);
    p0 = N(end);
  end
  N = setdiff(N,size(P,1)-(0:3));
  P = P(1:end-4,:);
end
