function [I,A,cen,vol] = inertia_tensor(V,F)
  % INERTIA_TENSOR Compute the inertia tensor and principal axes of a solid mesh
  % (V,F)
  % 
  % [I,A] = inertia_tensor(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of faces indices into F
  % Outputs:
  %   I  3x3 inertia tensor
  %   A  3x3 matrix with principal axes in rows
  %   cen  3d position of centroid
  %   vol  volume
  %
  % Adapted from RigidBodyParams by 
  % AUTHOR: Anton Semechko (a.semechko@gmail.com)
  % DATE  : December, 2014

  X=V;
  Tri=F;
  % Area weighted face normals (magnitude = 2 x triangle area)
  X1=X(Tri(:,1),:);
  X2=X(Tri(:,2),:);
  X3=X(Tri(:,3),:);
  FN=cross(X2-X1,X3-X1,2);
  
  % Zeroth order moment (same as volume enclosed by the mesh) ---------------
  C=(X1+X2+X3)/3;
  m000=sum(sum(FN.*C,2))/6;
  
  % First order moments (together they specify the centroid of the region 
  % enclosed by the mesh) ---------------------------------------------------
  x1=X1(:,1); y1=X1(:,2); z1=X1(:,3);
  x2=X2(:,1); y2=X2(:,2); z2=X2(:,3);
  x3=X3(:,1); y3=X3(:,2); z3=X3(:,3);
  
  x_2=((x1+x2).*(x2+x3) + x1.^2 + x3.^2)/12;
  y_2=((y1+y2).*(y2+y3) + y1.^2 + y3.^2)/12;
  z_2=((z1+z2).*(z2+z3) + z1.^2 + z3.^2)/12;
  
  xy=((x1+x2+x3).*(y1+y2+y3) + x1.*y1 + x2.*y2 + x3.*y3)/24;
  xz=((x1+x2+x3).*(z1+z2+z3) + x1.*z1 + x2.*z2 + x3.*z3)/24;
  yz=((y1+y2+y3).*(z1+z2+z3) + y1.*z1 + y2.*z2 + y3.*z3)/24;
  
  m100=sum(sum(FN.*[x_2 2*xy 2*xz],2))/6;
  m010=sum(sum(FN.*[2*xy y_2 2*yz],2))/6;
  m001=sum(sum(FN.*[2*xz 2*yz z_2],2))/6;
  
  % Second order moments (used to determine elements of the inertia tensor)
  % -------------------------------------------------------------------------
  x_3=((x1+x2+x3).*(x1.^2+x2.^2+x3.^2) + x1.*x2.*x3)/20; 
  y_3=((y1+y2+y3).*(y1.^2+y2.^2+y3.^2) + y1.*y2.*y3)/20; 
  z_3=((z1+z2+z3).*(z1.^2+z2.^2+z3.^2) + z1.*z2.*z3)/20; 
  
  x_2y=((3*y1+y2+y3).*x1.^2 + (y1+3*y2+y3).*x2.^2 + (y1+y2+3*y3).*x3.^2 + (2*y1+2*y2+y3).*x1.*x2 + (2*y1+y2+2*y3).*x1.*x3 + (y1+2*y2+2*y3).*x2.*x3)/60;
  x_2z=((3*z1+z2+z3).*x1.^2 + (z1+3*z2+z3).*x2.^2 + (z1+z2+3*z3).*x3.^2 + (2*z1+2*z2+z3).*x1.*x2 + (2*z1+z2+2*z3).*x1.*x3 + (z1+2*z2+2*z3).*x2.*x3)/60;
  
  y_2x=((3*x1+x2+x3).*y1.^2 + (x1+3*x2+x3).*y2.^2 + (x1+x2+3*x3).*y3.^2 + (2*x1+2*x2+x3).*y1.*y2 + (2*x1+x2+2*x3).*y1.*y3 + (x1+2*x2+2*x3).*y2.*y3)/60;
  y_2z=((3*z1+z2+z3).*y1.^2 + (z1+3*z2+z3).*y2.^2 + (z1+z2+3*z3).*y3.^2 + (2*z1+2*z2+z3).*y1.*y2 + (2*z1+z2+2*z3).*y1.*y3 + (z1+2*z2+2*z3).*y2.*y3)/60;
  
  z_2y=((3*y1+y2+y3).*z1.^2 + (y1+3*y2+y3).*z2.^2 + (y1+y2+3*y3).*z3.^2 + (2*y1+2*y2+y3).*z1.*z2 + (2*y1+y2+2*y3).*z1.*z3 + (y1+2*y2+2*y3).*z2.*z3)/60;
  z_2x=((3*x1+x2+x3).*z1.^2 + (x1+3*x2+x3).*z2.^2 + (x1+x2+3*x3).*z3.^2 + (2*x1+2*x2+x3).*z1.*z2 + (2*x1+x2+2*x3).*z1.*z3 + (x1+2*x2+2*x3).*z2.*z3)/60;
  
  xyz=((x1+x2+x3).*(y1+y2+y3).*(z1+z2+z3) - (y2.*z3+y3.*z2-4*y1.*z1).*x1/2 -(y1.*z3+y3.*z1-4*y2.*z2).*x2/2 - (y1.*z2+y2.*z1-4*y3.*z3).*x3/2)/60;
  
  m110=sum(sum(FN.*[x_2y y_2x 2*xyz],2))/6;
  m101=sum(sum(FN.*[x_2z 2*xyz z_2x],2))/6;
  m011=sum(sum(FN.*[2*xyz y_2z z_2y],2))/6;
  
  m200=sum(sum(FN.*[x_3 3*x_2y 3*x_2z],2))/9;
  m020=sum(sum(FN.*[3*y_2x y_3 3*y_2z],2))/9;
  m002=sum(sum(FN.*[3*z_2x 3*z_2y z_3],2))/9;
  
  % Inertia tensor ----------------------------------------------------------
  Ixx=m020+m002-(m010^2+m001^2)/m000;
  Iyy=m200+m002-(m100^2+m001^2)/m000;
  Izz=m200+m020-(m100^2+m010^2)/m000;
  Ixy=m110-m100*m010/m000;
  Ixz=m101-m100*m001/m000;
  Iyz=m011-m010*m001/m000;
  
  I=[Ixx -Ixy -Ixz; -Ixy Iyy -Iyz; -Ixz -Iyz Izz];
  
  % Output ------------------------------------------------------------------
  [sV,D]=svd(I);
  % Local frame of reference
  A=fliplr(sV);
  A(:,3)=cross(A(:,1),A(:,2));
  A = A';
  
  vol=m000;
  cen=[m100 m010 m001]/m000;
end
