function C = sample_interior(V,F,n)
  % SAMPLE_INTERIOR  Sample the interior of a solid bounded by the triangle
  % mesh (V,F)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   n  number of samples
  % Output:
  %   C  n by 3 list of samples
  %

  % SAMPLE_BOUNDING_BOX  Smample the interour of the bounding box of V
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   n  number of samples
  % Output:
  %   C  n by 3 list of samples
  %
  function C = sample_bounding_box(V,n)
    C = bsxfun(@plus, bsxfun(@times,rand(n,3),max(V)-min(V)),min(V));
  end

  out = 1:n;
  while ~isempty(out)
    C(out,:) = sample_bounding_box(V,numel(out));
    w_out = winding_number(V,F,C(out,:));
    out = out(abs(w_out)<0.5);
  end
end
