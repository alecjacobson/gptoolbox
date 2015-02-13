function S = snf(D)
  % SNF Smith normal form. This is the straight forward O(m*n*n)
  % implementation.
  %
  % S = snf(D)
  %
  % Inputs:
  %   D  m by n boundary operator matrix
  % Outputs:
  %   S  m by n Smith Normal Form matrix
  %

  function S = swap_rows(S,x,k)
    temp = S(x,:);
    S(x,:) = S(k,:);
    S(k,:) = temp;
  end
  function S = swap_cols(S,x,k)
    temp = S(:,x);
    S(:,x) = S(:,k);
    S(:,k) = temp;
  end

  S = logical(D);
  m = size(S,1);
  n = size(S,2);
  max_x = min([m n]);
  for x = 1:max_x
    [k,l] = find(S(x:end,x:end),1);
    if ~isempty(k)
      k = k+x-1;
      l = l+x-1;
      % exchange rows x and k 
      S = swap_rows(S,x,k);
      % exchange columns x and l 
      S = swap_cols(S,x,l);
      %% Not any faster, in fact, much slower on big input.
      %S(x+1:end,:) = mod(S(x+1:end,:)+bsxfun(@times,S(x+1:end,x),S(x,:)),2);
      %S(:,x+1:end) = mod(S(:,x+1:end)+bsxfun(@times,S(x,x+1:end),S(:,x)),2);
      for i = find(S(x+1:m,x))'
        i = i+x+1-1;
        % Slightly faster to use logicals and `xor` than doubles and `mod`
        S(i,:) = xor(S(i,:),S(x,:));
      end
      for j = find(S(x,x+1:n))
        j = j+x+1-1;
        S(:,j) = xor(S(:,j),S(:,x));
      end
    end
  end
end
