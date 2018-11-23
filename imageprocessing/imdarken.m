function C = imdarken(A,B)
  % IMDARKEN Blend A over B, replacing pixels in B if the cooresponding pixel in
  % A is darker.
  %
  % C = imdarken(A,B)
  %
  % Inputs:
  %   A  w * h * c image
  %   B  w * h * c image
  % Outputs:
  %   C  w * h * c image
  %  
  % Known issues:
  %   Does not incoorporate alpha channel like Photoshop would. 
  C = B;
  C(A<B) = A(A<B);
end
