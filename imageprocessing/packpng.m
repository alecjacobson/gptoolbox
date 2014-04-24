function rgba = packpng(G,filename,depth)
  % PACKPNG pack a grayscale double image into a depth bit RGBA PNG image
  %
  % rgba = packpng(G,filename,depth)
  %
  % Inputs:
  %   G  h by w image of doubles in [0,1]
  %   filename  path to png file
  %   depth  number of bits to use: 8 or 16
  % Outputs:
  %  rgba  RGBA png image
  %

  % clamp to [0,1]
  if any(G(:)<0)
    warning('Clamping values as small as %g to 0\n',min(G(:)));
    G(G<0) = 0;
  end
  if any(G(:)>1)
    warning('Clamping values as large as %g to 1\n',max(G(:)));
    G(G>1) = 1;
  end
  G_orig = G;

  % scale so 1 becomes (2^depth-1)/2^depth)
  G = G*(1 - 2^-depth);
  rgba = zeros(size(G,1),size(G,2),4);
  rgba_uint16 = uint16(zeros(size(G,1),size(G,2),4));
  for ii = 1:4
    rgba(:,:,ii) = floor(G*(2^depth))/(2^depth);
    rgba_uint16(:,:,ii) = uint16(rgba(:,:,ii)*2^16);
    G = (G - rgba(:,:,ii))*2^depth;
  end

  % imwrite is going to scale the values when we pass doubles, so we need to
  % preunscale them
  rgba = rgba/(1 - 2^-depth);
  imwrite(rgba(:,:,1:3),filename,'BitDepth',depth,'Alpha',rgba(:,:,4));
  [rgba2,~,alpha] = imread(filename);rgba2 = cat(3,rgba2,alpha);

end
