function [N,P,M,abs_M] = readIPI(filename)
  % READIPI read ipiSofts skeleton format
  % 
  % [N,P,M,abs_M] = readIPI(filename)
  %
  % Inputs:
  %   filename  path to ipiSoft file
  % Outputs:
  %   N  #B cell of bone names
  %   P  #B list of parent indices (0 means root) 
  %   M  4 by 4 by #B list of relative transform matrices
  %   abs_M 4 by 4 by #B list of absolute transform matrices
  %
  fp = fopen(filename,'r');
  line = eat_comments(fp,'#');
  fscanf(fp,'\n');
  while true
    ii = 1;
    % read id number
    line = line(ii:end);
    [id,count,e,ii] = sscanf(line,'%d ',1);
    if count ~= 1
      error('Bad format');
      break;
    end

    % read name
    line = line(ii:end);
    [name,count,e,ii] = sscanf(line,'%s ',1);
    if count ~= 1
      error('Bad format');
      break;
    end
    N{id} = name;

    % read parent
    line = line(ii:end);
    [p,count,e,ii] = sscanf(line,'%d ',1);
    if count ~= 1
      error('Bad format');
      break;
    end
    P(id) = p;

    % read parent
    line = line(ii:end);
    [m,count,e,ii] = sscanf(line,'%g ',16);
    if count ~= 16
      error('Bad format');
      break;
    end
    M(:,:,id) = reshape(m,4,4)';

    line = fgets(fp);
    if line == -1
      break;
    end
  end
  fclose(fp);

  next = P;
  abs_M = M;
  while true
    % loop over bones
    for ii = 1:numel(P)
      % if not root
      if next(ii) ~= 0
        abs_M(:,:,ii) = M(:,:,next(ii)) * abs_M(:,:,ii);
        next(ii) = P(next(ii));
      end
    end
    if ~any(next)
      break;
    end
  end

  assert(all(P<=numel(P)));
end
