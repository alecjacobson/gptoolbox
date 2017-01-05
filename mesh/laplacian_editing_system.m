function [A,rhs,pre_rhs,L_prime] = laplacian_editing_system(V,F,b,bc)
  % LAPLACIAN_EDITING_SYSTEM Build laplacian editing system matrix and right
  % hand side. Solution is obtained with sol = A\rhs or 
  % sol = A \ (pre_rhs * bc(:))
  %
  % [A,rhs,pre_rhs] = laplacian_editing_system(V,E,b,bc) 1D curve case
  % [A,rhs,pre_rhs] = laplacian_editing_system(V,F,b,bc) 2D surface case
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edge indices, [] means loop
  %   or
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   A  #V-#b by #V-b system matrix (**not** symmetric)
  %   rhs #V-#b by 1 right hand side
  %   pre_rhs #V-#b by #b pre right hand side multiplier matrix, so that
  %     rhs = pre_rhs*bc(:)
  %

  dim = size(V,2);
  assert(dim == 2);
  n = size(V,1);
  indices = 1:n;
  interior = indices(~ismember(indices,b));
  all = [interior(:);b(:)]';

  if(isempty(F) || size(F,2) ==2)
    % curve case
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Matlab script for 2D Laplacian Editing
    %
    % Please refer to "Laplacian Surface Editing" by Olga Sorkine,
    %   Daniel Cohen-Or, Yaron Lipman,  Marc Alexa, Christian Rössl and
    %   Hans-Peter Seidel,
    %   Eurographics/ACM SIGGRAPH Symposium on Geometry Processing 2004,
    %   pp. 179--188, ACM Press.
    %   http://www.cs.nyu.edu/~sorkine/ProjectPages/Editing/lse.html
    %
    % The script demonstrates how to build the Laplacian Editing matrix
    % when editing 2D curves. The static and handle vertices are hard-coded
    % for some example curves, provided with this script.
    % The script prompts the user to input new locations for the handle
    % vertices and displays the result in a new figure.
    %
    % Note: only the first stage of Laplacian editing is computed, namely
    % editing with implicit local similarity transformations. The second stage
    % (renormalizing the delta-coordinates to remove the introduced local
    % scaling) is not implemented here.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ONLY WORKS FOR IN-ORDER LOOPS
    if isempty(F)
      F = [1:n;2:n 1]';
    end

    % the Laplacian matrix (uniform weighting)
    A = adjacency_matrix(F);
    L = -0.5*sparse(diag(sum(A,2)) - A);

    delta = L*V;

    % we want to construct the matrix of the system for v-primes
    L_prime = [   L     sparse(n,n)   % the x-part
    	        sparse(n,n)    L    ]; % the y-part
    	      
    for i=1:n
      % the neighbors of i are i-1 and i+1
      ring = [find(A(i,:)) i];
      n_ring = numel(ring);
      ringV = V(ring,:)';
      ringV = [ringV; ones(1,n_ring)];

      % the coeff matrix for the system that solves for T
      %      s  a  t1
      % T = -a  s  t2
      %      0  0  1

      C = zeros(2*n_ring,4);
      % ... Fill C in
      for r=1:n_ring
        C(r,:) =                [ringV(1,r)       ringV(2,r)  ringV(3,r)      0  ];
        C(length(ring)+r,:) =   [ringV(2,r)  (-1)*ringV(1,r)       0  ringV(3,r) ];
      end;

      Cinv = pinv(C);
      s =  Cinv(1,:);
      a =  Cinv(2,:);

      delta_i = delta(i,:)';
      delta_ix = delta_i(1);
      delta_iy = delta_i(2);

      % T*delta gives us an array of coefficients
      Tdelta = [delta_ix*s      + delta_iy*a 
    	        delta_ix*(-1)*a + delta_iy*s];


      % updating the weights in Lx_prime, Ly_prime, Lw_prime
      L_prime(i,[ring (ring + n)]) = L_prime(i,[ring (ring + n)]) +...
                                                  (-1)*Tdelta(1,:);
      L_prime(i+n,[ring (ring + n)]) = L_prime(i+n,[ring (ring + n)]) +...
                                                    (-1)*Tdelta(2,:);
    end;

    % additional constraints - we need at least 3
    assert(numel(b) >= 3);

    A = ...
      L_prime([all (all+n)],[interior (interior+n)])' * ...
      L_prime([all (all+n)],[interior (interior+n)]);
    pre_rhs = ...
      - L_prime([all (all+n)],[interior (interior+n)])' * ...
      L_prime([all (all+n)],[b (b+n)]);
    rhs = pre_rhs * bc(:);
    %A = L_prime'*L_prime;
    %I = eye(size(A));
    %A([b (b+n)],:) = I([b (b+n)],:);
    %rhs = zeros(n*2,1);
    %rhs([b (b+n)]) = bc(:);
    %pre_rhs = I(:,[b (b+n)]);

  else
    % surface case
    assert(size(F,2) == 3);
    L = cotmatrix(V,F);
    delta = L * V;

    % we want to construct the matrix of the system for v-primes
    L_prime = [   L     zeros(n)   % the x-part
    	       zeros(n)    L    ]; % the y-part

    adj = adjacency_matrix(F);

    for i=1:n
      % the neighbors of i are i-1 and i+1
      ring = find(adj(i,:));
      ringV = V(ring,:)';
      ringV = [ringV; ones(1,length(ring))];

      n_ring = numel(ring);


      % the coeff matrix for the system that solves for T
      %      s  a  t1
      % T = -a  s  t2
      %      0  0  1

      C = zeros(2*n_ring,4);
      % ... Fill C in
      for r=1:n_ring
        C(r,:) =                [ringV(1,r)       ringV(2,r)  ringV(3,r)      0  ];
        C(n_ring+r,:) =   [ringV(2,r)  (-1)*ringV(1,r)       0  ringV(3,r) ];
      end;

      Cinv = pinv(C);
      s =  Cinv(1,:);
      a =  Cinv(2,:);

      delta_i = delta(i,:)';
      delta_ix = delta_i(1);
      delta_iy = delta_i(2);

      % T*delta gives us an array of coefficients
      Tdelta = [delta_ix*s      + delta_iy*a 
    	        delta_ix*(-1)*a + delta_iy*s];

    	    

      % updating the weights in Lx_prime, Ly_prime, Lw_prime
      L_prime(i,[ring (ring + n)]) = L_prime(i,[ring (ring + n)]) +...
                                                  (-1)*Tdelta(1,:);
      L_prime(i+n,[ring (ring + n)]) = L_prime(i+n,[ring (ring + n)]) +...
                                                    (-1)*Tdelta(2,:);
    end;

    % additional constraints - we need at least 3
    %assert(numel(b) >= 3);

    A = ...
      L_prime([all (all+n)],[interior (interior+n)])' * ...
      L_prime([all (all+n)],[interior (interior+n)]);
    pre_rhs = ...
      - L_prime([all (all+n)],[interior (interior+n)])' * ...
      L_prime([all (all+n)],[b (b+n)]);
    rhs = pre_rhs * bc(:);
    %A = L_prime'*L_prime;
    %I = eye(size(A));
    %A([b (b+n)],:) = I([b (b+n)],:);
    %rhs = zeros(n*2,1);
    %rhs([b (b+n)]) = bc(:);
    %pre_rhs = I(:,[b (b+n)]);
  end
end
