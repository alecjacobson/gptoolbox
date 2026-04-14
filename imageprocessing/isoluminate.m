function out = isoluminate(in,uniform)
  if nargin<2
    uniform = false;
  end
  size_in = size(in);
  rgb_in = reshape(in,[],3);
  lab_in = rgb2oklab(rgb_in);
  % find colors out such that:
  %
  %  min ‖ lab - lab_in ‖²
  %  such that lab(i,1) = lab(j,1) ∀ i,j
  %  and oklab2rgb(lab) ∈ [0,1]³
  %
  %  min ‖ [l*[1;1;1] ab] - lab_in ‖²
  %  min ‖ l - l_in ‖² + ‖ ab - ab_in ‖²
  %
  % Wait, is this boring? Suppose we know l★ then all l★ab lie on a plane
  % perpendicular to the line of grays. each l★ab's corresponding lab_in point
  % lies off this plane someplace, but the optimal l★ab will have the same
  % "argument"/"angle". So then ab★ will be as close to ab_in without going
  % outside the gamut. So either ab★=ab_in or ab★=max_arg(ab_in)
  %
  % So then it seems this is just a 1D search problem.
  % And it seems then mean(l) just minimizes it.

  if uniform
    min_E = inf;
    % brute-force search for best l
    for l = linspace(0,1)
      s = local(l,lab_in(:,2:3),0,1,true,0);
      lab = [repmat(l,size(lab_in,1),1) s.*lab_in(:,2:3)];
      E = sum((lab - lab_in).^2,'all');
      if E < min_E
        min_E = E;
        best_l = l;
        best_s = s;
      end
    end
    l = best_l;
    s = best_s;
  else
    l = mean(lab_in(:,1))
    s = local(l,lab_in(:,2:3),0,1,false,0);
  end
  lab = [repmat(l,size(lab_in,1),1) s.*lab_in(:,2:3)];
  rgb = oklab2rgb(lab);
  out = reshape(rgb,size_in);


  function s = local(l,ab_in,s_min,s_max,uniform,iter)
    if ~uniform && numel(s_min) == 1
      s_min = repmat(s_min,size(ab_in,1),1);
    end
    if ~uniform && numel(s_max) == 1
      s_max = repmat(s_max,size(ab_in,1),1);
    end
    if nargin<5;
      uniform = false;
    end
    if nargin<6;
      iter=0; 
    end
    if iter > 10
      s = s_min;
      return;
    end
    s = 0.5*(s_min+s_max);
    rgb = oklab2rgb([repmat(l,size(ab_in,1),1) s.*ab_in]);
    if uniform
      bad = any(rgb(:)<0 | rgb(:)>1);
    else 
      bad = any(rgb<0 | rgb>1,2);
    end
    s_max(bad)  =  s(bad);
    s_min(~bad) = s(~bad);
    s = local(l,ab_in,s_min,s_max,uniform,iter+1);
  end

end
