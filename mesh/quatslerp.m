function Qt = quatslerp(Q0,Q1,t)
  % QUATSLERP  Slerp between two quaterions
  %
  % Qt = quatslerp(Q0,Q1,t)
  % 
  % Inputs:
  %   Q0  #Q by 4 list of quaternions [w x y z] at t=0
  %   Q1  #Q by 4 list of quaternions at t=1
  %   t   #Q|1 list of times
  % Outputs:
  %   Qt  #Q by 4 list of quaternions at times t
  %

  if numel(t) == 1
    t = repmat(t,size(Q0,1),1);
  end

  % slerp pose
  d = sum(Q0.*Q1,2);
  absD = abs(d);
  theta = acos(absD);
  sinTheta = sin(theta);
  % https://github.com/alecjacobson/gptoolbox/issues/74
  s0 = sin( (1-t) .* theta)./sinTheta;
  s1 = sin(    t  .* theta)./sinTheta;
  s0(absD>=1) = 1-t(absD>=1);
  s1(absD>=1) =   t(absD>=1);
  s1(d<0) = -s1(d<0);
  Qt = bsxfun(@times,s0,Q0) + bsxfun(@times,s1,Q1);
  Qt = normalizerow(Qt);
end
