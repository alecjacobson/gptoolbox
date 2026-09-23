function M = isolines_stripe_map(M,val)
  if nargin < 2
    val = 0.9;
  end
  nm = size(M,1);
  des = repmat([1 1 1;val 1 1],ceil(nm/2),1);
  des = des(1:nm,:);
  M = max(min(oklab2rgb(des.*rgb2oklab(M)),1),0);

end
