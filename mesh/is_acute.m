function [IA,IO,e,A] = is_acute(V,F)
  % IS_ACUTE Determine if each triangle is a cute triangle.
  %
  % [IA,IO,A] = is_acute(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of face indices into V
  % Outputs:
  %   IA  #F list of flags revealing whether acute
  %   IO  #F list of flags revealing whether obtuse (not neessarily ~IO, but
  %     upto floating precision)
  %   e  #F by 3 list of edge lengths
  %   A  #F list of areas

  % computing upper bounds for each face in FA needs some auxilary quantities
  e = edge_lengths(V,F);
  % semiperimeter, area, circumradius, inradius
  s = sum(e,2)/2;
  A = doublearea(V,F);
  %A = sqrt(s.*(s-e(:,1)).*(s-e(:,2)).*(s-e(:,3)));
  R = e(:,1).*e(:,2).*e(:,3)./(4.0.*A);
  r = A./s;
  IA = s-r > 2*R;
  IO = s-r < 2*R;

end
