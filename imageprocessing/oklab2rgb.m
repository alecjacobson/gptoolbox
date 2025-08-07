function rgb = oklab2rgb(ok)
  % OKLAB2RGB Convert colors stored in OKLab color space to RGB color space
  %
  % rgb = oklab2rgb(ok)
  %
  % Inputs:
  %   ok  #h by #w by 3 image of colors or #ok by 3 list of colors
  % Outputs:
  %   rgb  converted ok colors (same size as ok)
  %
  % See also: rgb2lin, lin2rgb, oklab2lin , rgb2oklab, lin2oklab
  lin = oklab2lin(ok);
  if isreal(lin)
    rgb = lin2rgb(lin);
  else
    % Implement here with consideration of complex input so we can use complex
    % step finite differencing.
    x = lin;
    % From lin2rgb
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma = cast(1/2.4,'like',x);
    a     = cast(1.055,'like',x);
    b     = cast(-0.055,'like',x);
    c     = cast(12.92,'like',x);
    d     = cast(0.0031308,'like',x);
    y = zeros(size(x),'like',x);
    in_sign = -2 * (real(x) < 0) + 1;
    x = abs(real(x)) + imag(x)*1i;
    lin_range = (real(x) < d);
    gamma_range = ~lin_range;
    y(gamma_range) = a * exp(gamma .* log(x(gamma_range))) + b;
    y(lin_range) = c * x(lin_range);
    y = real(y) .* in_sign + imag(y)*1i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rgb = y;
  end
end

