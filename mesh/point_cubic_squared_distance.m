function [sqrD,S] = point_cubic_squared_distance(Q,C,tol)
  [P,T] = cubic_flat_eval(C,tol);
  E = [1:size(P,1)-1;2:size(P,1)]';
  [sqrD,I,C] = point_mesh_squared_distance(Q,P,E);
  B = barycentric_coordinates(C,P(E(I,1),:),P(E(I,2),:),'Project',true);
  S = T(E(I,1)).*B(:,1) + T(E(I,2)).*B(:,2);
end
