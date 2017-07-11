function [V,F] = schwarz_lantern(m,n)
  [X,Y,Z] = cylinder(ones(m,1),n);
  %[F,V] = surf2patch(X,Y,Z,'triangles');
  [Q,V] = surf2patch(X,Y,Z);
  R = axisangle2matrix([0 0 1],-2*pi/n*0.5);
  other = ((1:n+1)-1)*m+((2:2:m+1)');
  V(other,:) = V(other,:)*R;
  other = ((1:n)-1)*(m-1)+((1:2:m-1)');
  %tsurf(Q,V,'Tets',false,'FaceIndices',1,'CData',sparse(other,1,1,size(Q,1),1));
  %view(90,0);
  other = ismember(1:size(Q,1),other);
  F = [Q(~other,[1 2 3]);Q(~other,[1 3 4]);Q(other,[1 2 4]);Q(other,[4 2 3])];
  tsurf(F,V);
  %view(90,0);
end
