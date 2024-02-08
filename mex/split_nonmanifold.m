  % SPLIT_NONMANIFOLD Split a non-manifold (or non-orientable) mesh into a
  % manifold orientable mesh possibly with more connected components.
  %
  % Inputs:
  %   F  #F by 3 list of input tringle indices into some vertex list V
  %   %Optional:
  %   %  'V' followed by #V by dim vertex list
  % Outputs:
  %   SF  #F by 3 list of output tringle indices into V(SVI,:)
  %   SVI  #SV list of indices into V identifying vertex positions
  %   %SV  #SV by dim list of vertex positions so that SV = V(SVI,:)
  %
