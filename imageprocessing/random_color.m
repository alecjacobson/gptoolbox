function R = random_color(n,preset)
  % RANDOM_COLOR  generate a list of random colors
  %
  % R = random_color(size)
  % or
  % R = random_color(size,preset)
  %
  % Inputs:
  %   size  size of matrix of random colors to generate
  %   preset  string containing colormap name or 'Pastel'
  % Outputs:
  %   R  size by 3 list of RGB colors
  %

  
  if(~exist('n','var'))
    n = 1;
  end
  if(~exist('preset','var'))
    preset = '';
  end

  switch preset
      case ''
    R = rand([n,3]);
      case 'Pastel'
    attempt = 0;
    % probably should be between O(n) and O(log n)
    max_attempts = 10*prod(n);
    retry = true([prod(n),1]);
    R = zeros([prod(n),3]);
    while(attempt < max_attempts)
      num_retry = sum(retry);
      R(retry,:) = rand(num_retry,3);
      R(retry,:) = (291/255) * R(retry,:)./repmat(sum(R(retry,:),2),1,3);
      %retry(retry) = max(R(retry,:),[],2) > 0;
      % too bright
      next_retry = retry;
      next_retry(retry) = ...
        (max(R(retry,:),[],2) - min(R(retry,:),[],2)) > (221/255) ;
      % too grey
      next_retry(retry) = std(R(retry,:),0,2) < (68/255);
      if(~any(next_retry))
        break;
      end
      retry = next_retry;
      attempt = attempt+1;
    end
    R=reshape(R,[n 3]) +  34/255;
      case {'jet','hsv','hot','pink','flag','bone','gray','cool','copper'}
          f = str2func(preset);
          R = f(n*10);
          R = R(randperm(end),:);
          R = R(1:n,:);
      otherwise
          error('Unsupported preset');
  end

end
