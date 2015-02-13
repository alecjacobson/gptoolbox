function [p,theta,s,d,f] = screw(R1,t1,R2,t2)
  % SCREW Find a screw motion which interpolates two given rigid
  % transformations. In 2D, a screw motion is a rotation about a fixed point
  % (screw point) p. 
  %
  % Inputs:
  %   R1  dim by dim rotation matrix
  %   t1 1 by dim translation vector
  %   R2  dim by dim rotation matrix
  %   t2 1 by dim translation vector
  % Outputs:
  %   p  point on screw axis
  %   theta  angle of rotation around screw point/axis
  %   s  1 by dim screw axis (only in 3d)
  %   d  translation along screw axis
  %   f  1 by dim fudge factor translation (if change in rotation is 0)
  %
  % Example:
  %   % get a polygon
  %   P = getPosition(impoly);
  %   E = [1:size(P,1);2:size(P,1) 1]';
  %   % Helper functions
  %   rotation = @(theta) [cos(theta) -sin(theta);sin(theta) cos(theta)];
  %   around = @(P,R,p) bsxfun(@plus,bsxfun(@minus,P,p)*R,p);
  %   rigid = @(P,R,t) bsxfun(@plus,P*R,t);
  %   R1 = rotation(pi/3);
  %   t1 = [0 1];
  %   R2 = rotation(2*pi/3);
  %   t2 = [1 0];
  %   [p,theta,~,~,f] = screw(R1,t1,R2,t2);
  %   plot_edges(rigid(P,R1,t1),E,'g','LineWidth',4);
  %   hold on;
  %   plot_edges(rigid(P,R2,t2),E,'r','LineWidth',4);
  %   scatter(p(1),p(2),'ok','LineWidth',4);
  %   for t = linspace(0,1,30)
  %     Q = bsxfun(@plus,around(rigid(P,R1,t1),rotation(t*theta),p),t*f);
  %     plot_edges(Q,E); 
  %   end
  %   hold off;
  %   axis equal;
  %
  % Example:
  %   R1 = axisangle2matrix([0 0 1],pi/3);
  %   R2 = axisangle2matrix([0 0 1],2*pi/3);
  %   t1 = 0*[1 2 -1];
  %   t2 = 0*[-2 1 -1];
  %   rigid = @(P,R,t) bsxfun(@plus,P*R,t);
  %   around = @(P,theta,p,s) ...
  %     bsxfun(@plus,bsxfun(@minus,P,p)*axisangle2matrix(s,theta),p);
  %   tsurf(F,rigid(V,R1,t1),'FaceColor','g');
  %   [p,theta,s,d,f] = screw(R1,t1,R2,t2);
  %   hold on;
  %   tsurf(F,rigid(V,R2,t2),'FaceColor','r');
  %   for t = linspace(0.1,0.9,8);
  %     U = bsxfun(@plus,around(rigid(V,R1,t1),theta*t,p,s),t*d*s+t*f);
  %     tsurf(F,U,'CData',(1:size(F,1))');
  %   end
  %   quiver3(p(1),p(2),p(3),s(1),s(2),s(3),'ok','LineWidth',4);
  %   hold off;
  %   axis equal;
  %  

  dim = size(R1,1);
  assert(all(size(R1)==dim));
  assert(all(size(R2)==dim));
  assert(numel(t1)==dim);
  assert(numel(t2)==dim);
  % row vectors
  t1 = reshape(t1,[],dim);
  t2 = reshape(t2,[],dim);

  switch dim
  case 2
    % Unused in 2d, set as if in 3d
    s = [0 0 1];
    d = 0;
    if norm(R2-R1,'fro')<eps
      f = t2-t1;
      p = [0 0];
      theta = 0;
    else
      % invert basis according to x*R1+t
      t = (t2-t1)*R1';
      R = R2*R1';
      % Fixed point of transformation
      % pR+t = p
      % p(R-I) = -t
      I = eye(dim);
      p = (-t)/(R-I);
      p = p*R1+t1;
      % Change basis back according to x*R1+t
      theta = atan2(R2(2),R2(1)) - atan2(R1(2),R1(1));
      if theta > pi
        theta = -(2*pi - theta);
      end
      f = [0 0];
    end
  case 3
    % "Collision Prediction for Polyhedra under Screw Motions" [Kim & Rossignac
    % 03]
    dR = R2-R1;
    if sum(sum(dR.^2,2)<eps)>=2
      f = t2-t1;
      p = [0 0 0];
      s = [0 0 1];
      theta = 0;
      d = 0;
    else
      s_tilde = ...
        cross(dR(1,:),dR(2,:),2) + ...
        cross(dR(2,:),dR(3,:),2) + ...
        cross(dR(3,:),dR(1,:),2);
      s = normalizerow(s_tilde);
      % Should be able to rewrite this using atan2
        sqrt(sum((dR(1,:)).^2,2))
        sqrt(sum(cross(s,R1(1,:),2).^2,2))
      theta = 2 * asin( ...
        0.5*sqrt(sum((dR(1,:)).^2,2))./ ...
        sqrt(sum(cross(s,R1(1,:),2).^2,2)));
      % Just consider the point o as the origin
      o1 = bsxfun(@plus,[0 0 0]*R1,t1);
      o2 = bsxfun(@plus,[0 0 0]*R2,t2);
      d = sum((o2-o1).*s,2);
      % Where o' where did this sign flip come from?
      p = 0.5 * (o2 + o1 + cross(s,o2-o1,2)./tan(-theta/2));
      f = [0 0 0];
    end
  end

end
