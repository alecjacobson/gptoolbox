function D = div_intrinsic(l,F)
  % DIV_INTRINSIC Construct an intrinsic divergence operator.
  %
  % Inputs:
  %  l  #F by 3 list of edge lengths
  %  F  #F by 3 list of triangle indices into some vertex list V
  % Outputs:
  %  G  #V by #F*2 divergence matrix.
  %
  % See also: grad_intrinsic
  %

  G = grad_intrinsic(l,F);
  TA = repdiag(sparse(diag(doublearea_intrinsic(l))),2);
  D = -0.25*G'*TA;

end

