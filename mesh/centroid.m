function [C,vol] = centroid(V,F,varargin)
  % CENTROID Compute the centroid of a closed polyhedron boudned by (V,F)
  %
  % C = centroid(V,F)
  % [C,vol] = centroid(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %   Optional:
  %     'Robust' followed by whether to use more robust but costlier method for
  %       nearly closed input. {false}
  % Outputs:
  %   C  3-vector of centroid location
  %   vol  total volume of polyhedron
  %

  % default values
  robust = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Robust'}, {'robust'});
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

  cross2 = @(a,b,c) ...
    [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
     a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];

  if robust
    assert(size(V,2) == 3,'Only 3d supported');
    [TV,TT] = cdt(V,F);
    BC = barycenter(TV,TT);
    w = winding_number(V,F,barycenter(TV,TT))/(4*pi);
    TT = TT(abs(w)>0.5,:);
    BC = BC(abs(w)>0.5,:);
    vol = volume(TV,TT);
    C = sum(bsxfun(@times,vol,BC))/sum(vol);
    vol = sum(vol);
  else
    % "Calculating the volume and centroid of a polyhedron in 3d" [Nuernberg 2013]
    % http://www2.imperial.ac.uk/~rn/centroid.pdf
    switch size(V,2)
    case 2
      % https://en.wikipedia.org/wiki/Centroid#Centroid_of_polygon
      % Rename corners
      A = V(F(:,1),:);
      B = V(F(:,2),:);
      D = A(:,1).*B(:,2) - B(:,1).*A(:,2);
      vol = 0.5*sum(D);
      C = (1/(6*vol)).*sum((A+B).*D);
    case 3
      % Rename corners
      A = V(F(:,1),:);
      B = V(F(:,2),:);
      C = V(F(:,3),:);
      % Needs to be **unnormalized** normals
      N = cross2(B-A,C-A);
      % total volume via divergence theorem: ∫ 1
      vol = sum(sum(A.*N))/6;
      % centroid via divergence theorem and midpoint quadrature: ∫ x
      C = 1/(2*vol)*(1/24* sum(N.*((A+B).^2 + (B+C).^2 + (C+A).^2)));
    end
  end

end
