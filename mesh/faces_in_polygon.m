function [xin,iin,cin] = faces_in_polygon(V,F,X,Y,Vin)
  % FACES_IN_POLYGON test whether faces of mesh are inside a given polygon
  %  
  % [xin,iin,cin] = faces_in_polygon(V,F,X,Y)
  % [xin,iin,cin] = faces_in_polygon(V,F,X,Y)
  % [xin,iin,cin] = faces_in_polygon(V,F,X,Y,Vin)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of face indices
  %   X  #P by 1 list of polygon x-coordinates
  %   Y  #P by 1 list of polygon y-coordinates
  %   Optional
  %     Vin  #V list of flags overriding whether each position in V is defined
  %       as inside polygon formed by (X,Y), default is to use inpolygon
  % Outputs:
  %   xin #F list of flags revealing whether query face indices are all in (X,Y)
  %   iin #F list of flags revealing whether query faces have at least one
  %     index in (X,Y)
  %   cin #F list of flags revealing whether query face barycenters are in (X,Y)
  %
  % See in_mesh, inpolygon
  %

  % NaNs are not recognized as they are in inpolygon
  assert(~any(isnan(X)));
  assert(~any(isnan(Y)));

  dim = size(V,2);
  % only works in 2D
  assert(dim == 2);

  if ~exist('exclusive','var')
    exclusive = true;
  end

  if ~exist('Vin','var')
    % first determine mesh points strictly in or on polygon
    Vin = inpolygon(V(:,1),V(:,2),X,Y);
  end

  % find faces that are inside because all vertices are inside
  [~,xin] = limit_faces(F,Vin,true);
  [~,iin] = limit_faces(F,Vin,false);

  % Note: if you can assume that if a triangle's corners are all inside (X,Y)
  % then also it's centroid is inside, then you could speed this up by only
  % testing triangles not in xin

  % barycenters for all triangles
  B = barycenter(V,F);

  % only do this if asked to
  if nargout >= 3
    % determine if each triangle barycenter is in polygon (X,Y)
    cin = inpolygon(B(:,1),B(:,2),X,Y);
  end



end
