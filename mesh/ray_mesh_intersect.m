function [flag, t, lambda] = ray_mesh_intersect(o, d, V, F);
% RAY_MESH_INTERSECT  Ray/mesh intersection using the algorithm proposed by
% Möller and Trumbore (1997).
%
% [flag, t, lambda] = ray_mesh_intersect(o, d, V, F);
%
% Input:
%    o  3D vector ray origin.
%    d  3D vector ray direction.
%    V  #V by 3 list of vertex positions
%    F  #F by 3 list of triangle indices
% Output:
%    flag  #F list of bools: (false) Reject, (true) Intersect.
%    t  #F list of distances from the ray origin.
%    lambda: #F by 3 list of barycentric coordinate of hits on each triangle
% Modified from code by: 
%    Jesus Mena

    epsilon = eps;%0.00001;

    assert(size(V,2) == 3);
    assert(size(F,2) == 3);
    % number of triangles
    m = size(F,1);
    assert(numel(d) == 3);
    % make direction a row vector
    d = reshape(d,1,3);
    %d = normalizerow(d);
    d = d ./ sqrt(sum(d.^2,2));
    assert(numel(o) == 3);
    % make origin a row vector
    o = reshape(o,1,3);

    p1 = V(F(:,1),:);
    p2 = V(F(:,2),:);
    p3 = V(F(:,3),:);

    %edge vector from corner 3 to 1
    e31 = p1-p3;
    %edge vector from corner 3 to 2
    e32 = p2-p3;

    % make copy of ray for each triangle
    %D = repmat(d,m,1);
    %O = repmat(o,m,1);

    %q  = cross(D,e32,2);
    %q  = cross2(D,e32);
    q  = crossrow(e32,-d);
    %a  = dot(e31,q,2); % determinant of the matrix M
    a  = dot2(e31,q); % determinant of the matrix M

    % I THINK THIS CAN BE SPED UP A LOT BY USING SPARSE VECTORS THROUGHOUT AND
    % NOT COMPUTING t,lambda FOR MISSES

    % assume all are hits
    flag = true(m,1);
    lambda = zeros(m,3);

    flag(a>-epsilon & a<epsilon) = 0;
    
    f = 1./a;
    %s = O-p3;
    %s = -(-O+p3);
    s = -plusrow(p3,-o);
    %lambda(:,1) = f.*dot(s,q,2);
    lambda(:,1) = f.*dot2(s,q);
    

    % the intersection is outside of the triangle
    flag(lambda(:,1)<0.0)  = 0;
    
    %r = cross(s,e31,2);
    r = cross2(s,e31);
    %lambda(:,2) = f.*dot(D,r,2);
    %lambda(:,2) = f.*dot2(D,r);
    lambda(:,2) = f.*dotrow(r,d);
    
    % the intersection is outside of the triangle
    flag(lambda(:,2)<0.0 | (lambda(:,1)+lambda(:,2))>1.0) = 0;

    %t = f.*dot(e32,r,2); % verified! 
    t = f.*dot2(e32,r); % verified! 
    flag(t < 0) = 0;

  function r = dotrow(a,b)
    % Computes dot product of rows in a against b.
    % optimizes r = dot(a,repmat(b,size(a,1),1),2)
    r = a(:,1).*b(1,1) + a(:,2).*b(1,2) + a(:,3).*b(1,3);
  end

  function r = crossrow(a,b)
    % Computes cross product of rows in a against b.
    % optimizes r = cross(a,repmat(b,size(a,1),1),2)
    r =[a(:,2).*b(1,3)-a(:,3).*b(1,2), ...
        a(:,3).*b(1,1)-a(:,1).*b(1,3), ...
        a(:,1).*b(1,2)-a(:,2).*b(1,1)];
  end

  function r = plusrow(a,b)
    % Computes sum rows in a with b.
    % optimizes r = a + repmat(b,size(a,1),1)
    r = [a(:,1)+b(1,1) a(:,2)+b(1,2) a(:,3)+b(1,3)];
  end

  function r = dot2(a,b)
    % Optimizes r = dot(a,b,2), that is it computes dot products per row
    % Faster than dot if I know that I'm calling it correctly
    r = sum(a.*b,2);
  end
  function r = cross2(a,b)
    % Optimizes r = cross(a,b,2), that is it computes cross products per row
    % Faster than cross if I know that I'm calling it correctly
    r =[a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
        a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
        a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  end
end

