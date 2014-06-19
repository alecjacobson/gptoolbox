function [theta,cos_theta] = dihedral_angles(V,T,varargin)
  % DIHEDRAL_ANGLES Compute dihedral angles for all tets of a given tet mesh
  % (V,T)
  %
  % theta = dihedral_angles(V,T)
  % theta = dihedral_angles(V,T,'ParameterName',parameter_value,...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   T  #V by 4 list of tet indices
  %   Optional:
  %     'SideLengths' followed by #T by 6 list of tet edges lenghts: 
  %       41 42 43 23 31 12
  %     'FaceAreas' followed by #T by 4 list of tet face areas
  % Outputs:
  %   theta  #T by 6 list of dihedral angles (in radians)
  %   cos_theta  #T by 6 list of cosine of dihedral angles (in radians)
  %

  l = [];
  s = [];
  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'SideLengths'
      assert((v+1)<=numel(varargin));
      v = v+1;
      l = varargin{v};
    case 'FaceAreas'
      assert((v+1)<=numel(varargin));
      v = v+1;
      s = varargin{v};
    otherwise
      error(['Unsupported parameter: ' varargin{v}]);
    end
    v = v+1;
  end
  if isempty(l)
    % lengths of edges opposite *face* pairs: 23 31 12 41 42 43
    l = [ ...
      sqrt(sum((V(T(:,4),:)-V(T(:,1),:)).^2,2)) ...
      sqrt(sum((V(T(:,4),:)-V(T(:,2),:)).^2,2)) ...
      sqrt(sum((V(T(:,4),:)-V(T(:,3),:)).^2,2)) ...
      sqrt(sum((V(T(:,2),:)-V(T(:,3),:)).^2,2)) ...
      sqrt(sum((V(T(:,3),:)-V(T(:,1),:)).^2,2)) ...
      sqrt(sum((V(T(:,1),:)-V(T(:,2),:)).^2,2)) ...
    ];
  end
  % (unsigned) face Areas (opposite vertices: 1 2 3 4)
  if isempty(s)
    s = 0.5*[ ...
      doublearea_intrinsic(l(:,[2 3 4])) ...
      doublearea_intrinsic(l(:,[1 3 5])) ...
      doublearea_intrinsic(l(:,[1 2 6])) ...
      doublearea_intrinsic(l(:,[4 5 6]))];
  end

  % Law of cosines
  % http://math.stackexchange.com/a/49340/35376
  H_sqr = (1/16) * ...
    (4*l(:,[4 5 6 1 2 3]).^2.* l(:,[1 2 3 4 5 6]).^2 - ...
    ((l(:, [2 3 4 5 6 1]).^2 + l(:,[5 6 1 2 3 4]).^2) - ...
     (l(:, [3 4 5 6 1 2]).^2 + l(:,[6 1 2 3 4 5]).^2)).^2);
  cos_theta= (H_sqr - s(:,[2 3 1 4 4 4]).^2 - ...
                      s(:,[3 1 2 1 2 3]).^2)./ ...
                  (-2*s(:,[2 3 1 4 4 4]).* ...
                      s(:,[3 1 2 1 2 3]));
  % If this line ever becomes a bottleneck (unlikely) we should this function
  % as a wrapper to cos_dihedral_angles(V,T) or better as a parameter 'CosOnly',true
  theta = acos(cos_theta);

  % Q: Could also use atan2 formulation if you're after the dihedral angles (not
  % the cosine of dihedral angles)?
  % A: Yup, http://en.wikipedia.org/wiki/Dihedral_angle#Alternative_definitions
end
