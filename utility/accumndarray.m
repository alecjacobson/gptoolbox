function H = accumndarray(varargin)
  % H = accumndarray(I1,I2,...,IN,data,sz)
  %
  % Subscripted and broadcasted version of accumarray.
  %
  % I1, I2, ..., IN, and data can be nd-arrays of matching or "broadcastable"
  % sizes. Where "broadcastable" means that sizes that don't match are 1. 
  % 
  % If sizes match then this is equivalent to:
  % 
  % H = accumarray([I(:) J(:) K(:) ...],data(:),sz);
  %
  % N should be equal to numel(sz)
  %
  %
  vec = @(X) reshape(X,[],1);
  sz = varargin{end};
  subs = varargin(1:end-2);

  dsub = [];
  for si = 1:numel(subs)
    di = size(subs{si});
    dsub(end+1:numel(di)) = 1;
    dsub(1:numel(di)) = max(dsub(1:numel(di)),di);
  end
  data = varargin{end-1};
  for si = 1:numel(subs)+1
    if si>numel(subs)
      di = size(data);
    else
      di = size(subs{si});
    end
    di(end+1:numel(dsub)) = 1;
    assert(all(di == 1 | di == dsub));
    ri = dsub;
    ri(di==dsub) = 1;
    if si>numel(subs)
      data = vec(repmat(data,ri));
    else
      subs{si} = vec(repmat(subs{si},ri));
    end
  end
  H = accumarray([subs{:}],data,sz);
end

