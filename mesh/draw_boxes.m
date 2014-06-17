function p = draw_boxes(BN,BX,varargin)
  % DRAW_BOXES  Draw a list of boxes given minimum corners and maximum corners
  % 
  % p = draw_boxes(BN,BX)
  % 
  % Inputs:
  %   BN  #B by dim list of minimum corners
  %   BX  #B by dim list of maximum corners
  % Outputs:
  %   p  plot handle
  %
  dim = size(BN,2);
  m = size(BN,1);
  switch dim
  case 3
    V = [];
    for x = 0:1
      for y = 0:1
        for z = 0:1
          V = [V;BN + bsxfun(@times,BX-BN,[x y z])];
        end
      end
    end
    F = [3 4 2 1 
         1 2 6 5 
         8 6 2 4 
         4 8 7 3 
         5 7 3 1
         5 7 8 6
         ];
    F = reshape(bsxfun(@plus,(F(:)-1)*m,1:m)',6*m,4);
    t = trisurf(F,V(:,1),V(:,2),V(:,3),varargin{:});
  end
end
