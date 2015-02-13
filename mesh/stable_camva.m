function stable_camva(a)
  % STABLE_CAMVA Act like camva without changing zoom (aka dolly zoom)
  %
  % Inputs:
  %   a  angle input to camva
  %
  % See also: camva
  old_a = camva;
  old_t = camtarget;
  old_p = campos;
  camva(a);
  s = tan(a/2*pi/180)/tan(old_a/2*pi/180);
  p = (old_p-old_t)/s+old_t;
  campos(p);
end
