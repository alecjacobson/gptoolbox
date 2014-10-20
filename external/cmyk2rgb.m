function rgb = cmyk2rgb(cmyk)
  rgb = applycform(cmyk,makecform('cmyk2srgb'));
end
