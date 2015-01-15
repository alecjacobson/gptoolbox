function C = cotangent(V,F,varargin)
  % COTANGENT compute the cotangents of each angle in mesh (V,F), more details
  % can be found in Section 1.1 of "Algorithms and Interfaces for Real-Time
  % Deformation of 2D and 3D shapes" [Jacobson 2013]
  % 
  % C = cotangent(V,F)
  % C = cotangent(V,F,'ParameterName',parameter_value,...)
  %
  % Known bugs:
  %   This seems to return 0.5*C and for tets already multiplies by
  %   edge-lengths
  %
  % Inputs:
  %   V  #V by dim list of rest domain positions
  %   F  #F by {3|4} list of {triangle|tetrahedra} indices into V
  %   Optional (3-manifolds only):
  %     'SideLengths' followed by #F by 3 list of edge lengths corresponding
  %       to: 23 31 12. In this case V is ignored.
  %       or
  %       followed by #T by 6 list of tet edges lengths:
  %       41 42 43 23 31 12
  %     'FaceAreas' followed by #T by 4 list of tet face areas
  % Outputs:
  %   C  #F by {3|6} list of cotangents corresponding
  %     angles for triangles, columns correspond to edges 23,31,12
  %     dihedral angles *times opposite edge length* over 6 for tets, 
  %       WRONG: columns correspond to edges 23,31,12,41,42,43 
  %       RIGHT: columns correspond to *faces* 23,31,12,41,42,43
  %
  % See also: cotmatrix
  %
  % Copyright 2013, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  switch size(F,2)
  case 3

    % default values
    l = [];
    % Map of parameter names to variable names
    params_to_variables = containers.Map( ...
      {'SideLengths'}, {'l'});
    v = 1;
    while v <= numel(varargin)
      param_name = varargin{v};
      if isKey(params_to_variables,param_name)
        assert(v+1<=numel(varargin));
        v = v+1;
        % Trick: use feval on anonymous function to use assignin to this workspace 
        feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
      else
        error('Unsupported parameter: %s',varargin{v});
      end
      v=v+1;
    end
    assert(numel(varargin) == 0);

    % triangles
    if isempty(l)
      % edge lengths numbered same as opposite vertices
      l = [ ...
        sqrt(sum((V(F(:,2),:)-V(F(:,3),:)).^2,2)) ...
        sqrt(sum((V(F(:,3),:)-V(F(:,1),:)).^2,2)) ...
        sqrt(sum((V(F(:,1),:)-V(F(:,2),:)).^2,2)) ...
        ];
    end
    i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);
    l1 = l(:,1); l2 = l(:,2); l3 = l(:,3);
    % semiperimeters
    s = (l1 + l2 + l3)*0.5;
    % Heron's formula for area
    dblA = 2*sqrt( s.*(s-l1).*(s-l2).*(s-l3));
    % cotangents and diagonal entries for element matrices
    % correctly divided by 4 (alec 2010)
    C = [ ...
      (l2.^2 + l3.^2 -l1.^2)./dblA/4 ...
      (l1.^2 + l3.^2 -l2.^2)./dblA/4 ...
      (l1.^2 + l2.^2 -l3.^2)./dblA/4 ...
      ];
  case 4
    % Dealing with tets
    T = F;

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

    [~,cos_theta] = dihedral_angles([],[],'SideLengths',l,'FaceAreas',s);

    % volume
    vol = volume_intrinsic(l);

    %% Law of cosines
    %% http://math.stackexchange.com/a/49340/35376
    %H_sqr = (1/16) * ...
    %  (4*l(:,[4 5 6 1 2 3]).^2.* l(:,[1 2 3 4 5 6]).^2 - ...
    %  ((l(:, [2 3 4 5 6 1]).^2 + l(:,[5 6 1 2 3 4]).^2) - ...
    %  (l(:, [3 4 5 6 1 2]).^2 + l(:,[6 1 2 3 4 5]).^2)).^2);
    %cos_theta= (H_sqr - s(:,[2 3 1 4 4 4]).^2 - ...
    %                    s(:,[3 1 2 1 2 3]).^2)./ ...
    %                (-2*s(:,[2 3 1 4 4 4]).* ...
    %                    s(:,[3 1 2 1 2 3]));
    %% To retrieve dihedral angles stop here...
    %theta = acos(cos_theta);
    %theta/pi*180

    % Law of sines
    % http://mathworld.wolfram.com/Tetrahedron.html
    %sin_theta = bsxfun(@rdivide,vol,(2./(3*l)) .* s(:,[1 2 3 2 3 1]) .* s(:,[4 4 4 3 1 2]));
    sin_theta = bsxfun(@rdivide,vol,(2./(3*l)) .* s(:,[2 3 1 4 4 4]) .* s(:,[3 1 2 1 2 3]));
    %% Using sin for dihedral angles gets into trouble with signs
    %theta = asin(sin_theta);
    %theta/pi*180

    % http://arxiv.org/pdf/1208.0354.pdf Page 18
    C = 1/6 * l .* cos_theta ./ sin_theta;
    %% One should probably never use acot to retrieve signed angles
    %theta = acot(C);
    %theta/pi*180
  otherwise
    error('Unsupported simplex type');
  end

end
