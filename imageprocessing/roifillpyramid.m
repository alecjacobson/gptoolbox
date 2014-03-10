function J = roifillpyramid(I,BW)
  % ROIFILLPYRAMID reimplement matlab's image processing toolbox's roifill
  % function using multiresolution image pyramid and iterative solve for
  % laplace equation.
  %
  % J = roifillpyramid(I,BW)
  %
  % Inputs:
  %   I  h by w by c source image with hole
  %   BW  h by w binary image specifying where to fill
  % Outputs:
  %   J  h by w by c Filled image
  %

  assert(size(I,1) == size(BW,1));
  assert(size(I,2) == size(BW,2));
  % width
  w = size(I,2);
  % height
  h = size(I,1);
  % number of channels
  nc = size(I,3);

  nw = 2^ceil(log2(max([w h])))+1;
  nh = nw;
  NI = zeros( [ nw nh nc ]);
  NI(1:h,1:w,:) = I;
  NBW = false( [ nw nh nc ]);
  NBW(1:h,1:w,:) = repmat(BW,[1 1 nc]);


  % number of levels in pyramid
  nl = floor(log2(min([nw h])));
  % pyramid image of original image
  PI = {};
  PBW = {};
  PI{1} = NI;
  % pyramid image of original mask
  PBW{1} = NBW;
  % reduce until at smallest level
  for( level = 2:nl)
    PI{level} = impyramid(PI{level-1},'reduce');
    PBW{level} = impyramid(PBW{level-1},'reduce');
  end
  f = fspecial('gaussian',[3 3]);
  % expand until back at largest level
  for( level = ((nl:-1:2)))

    up = impyramid(PI{level},'expand');
    old = PI{level-1};
    PI{level-1} = ...
      (~PBW{level-1}).*PI{level-1} + ...
      ( PBW{level-1}).*up;

    if(any(PBW{level-1}(:)))
      for( c = 1:nc)
        diff = 1;
        min_diff = 1e-5;
        iter = 0;
        max_iterations = 100;
        while(diff > min_diff && iter < max_iterations)
          old = PI{level-1}(:,:,c);
          PI{level-1}(:,:,c) = roifilt2(f,PI{level-1}(:,:,c),PBW{level-1}(:,:,c));
          diff = sqrt(sum(sum((old-PI{level-1}(:,:,c)).^2)))/size(PBW{level-1},1);
          iter = iter + 1;
        end
        if(iter >= max_iterations)
          fprintf('Max iterations reached...\n');
        end
        if(diff < min_diff)
          fprintf('Minimal change reached...\n');
        end
      end
    end

  end

  J = PI{1}(1:h,1:w,:);
end
