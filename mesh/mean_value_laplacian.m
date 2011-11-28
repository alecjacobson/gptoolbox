function L = mean_value_laplacian(V,F)
  % MEAN_VALUE_LAPLACIAN Discrete laplacian using mean value weights. Notice
  % that this laplacian will not be symmetric or positive semi-definite.
  %
  % L = mean_value_laplacian(V,F)
  %
  % Inputs:
  %  V: #V x dim matrix of vertex coordinates
  %  F: #F by 3, list of indices of triangle corners
  % Outputs:
  %   L  sparse nvert x nvert matrix of mean value weights
  %
  % See also: cotmatrix, wachspress_laplacian

  % number of vertices
  n = size(V,1);

  % edge lengths numbered same as opposite vertices
  l = [ ...
    sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
    sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
    sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
    ];

  % should change code below, so we don't need this transpose
  if(size(F,1) == 3 && size(F,2) ~= 3)
    warning('F seems to be 3 by #F, it should be #F by 3');
  end

  % 
  %     vj+1
  %     /  \
  %    /    \
  %   /βij   \
  % vi--------vj
  %   \αij   /
  %    \    /
  %     \  /
  %     vj-1
  % 
  %
  % Lij = tan(αij/2) + tan(βij/2) )/|vi -vj|^2

  % http://en.wikipedia.org/wiki/Tangent_half-angle_formula
  % tan(θ/2) = sin θ/(1+cos θ)

  % renaming indices of vertices of triangles for convenience
  i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); 
  l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);

  % semiperimeters
  s = (l1 + l2 + l3)*0.5;
  % Heron's formula for area
  dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
  halftan12 = 2*dblA ./ (2*l1.*l2+l1.^2+l2.^2-l3.^2);
  halftan23 = 2*dblA ./ (2*l2.*l3+l2.^2+l3.^2-l1.^2);
  halftan31 = 2*dblA ./ (2*l3.*l1+l3.^2+l1.^2-l2.^2);

  % throw each angle at each edge entry
  Li = [ ...
    i1;i1;...
    i2;i2;...
    i3;i3;...
    ];
  Lj = [ ...
    i2;i3;...
    i3;i1;...
    i1;i2;...
    ];
  % divide by edge lengths
  Lv = [ ...
    halftan23./(l3);halftan23./(l2);...
    halftan31./(l1);halftan31./(l3);...
    halftan12./(l2);halftan12./(l1);...
    ];
  %Lv = [ ...
  %  halftan23;halftan23;...
  %  halftan31;halftan31;...
  %  halftan12;halftan12;...
  %  ];
  L = sparse(Li,Lj,Lv,n,n);
  %% normalize weights
  %L = bsxfun(@rdivide,L,sum(L,2));
  % set diagonal entries
  L = L - diag(sum(L,2));
end
