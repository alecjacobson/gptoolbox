function cmyk = rgb2cmyk(rgb)
  cmyk = applycform(rgb,makecform('srgb2cmyk'));
end

