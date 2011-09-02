function imshow_before_and_after(before_R,before_G,before_B,after_R,after_G,after_B,x,y,w,h)
  % show a clip of two images x,y,w,h side by side in a window
  if(nargin<7)
    x = 1; y = 1; w = size(before_R,2); h = size(before_R,1);
  end

  imshow_rgb( ...
    [before_R(y:y+h,x:x+h), after_R(y:y+h,x:x+w)], ...
    [before_G(y:y+h,x:x+h), after_G(y:y+h,x:x+w)], ...
    [before_B(y:y+h,x:x+h), after_B(y:y+h,x:x+w)]);
end
