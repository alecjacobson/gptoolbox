function [P,C] = parse_path(dstr)
  % [P,C] = parse_path(dstr)
  %
  % Inputs:
  %   dstr  string of the value of the 'd' attribute
  % Outputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %
  function [x,dstr] = parse_x(dstr)
    [x,count,~,pos] = sscanf(dstr,'%g',1);
    if count~= 1
      xy = [];
      return;
    end
    dstr = dstr(pos:end);
    dstr = strip(strip(dstr,'left',','),'left',' ');
  end
  function [xy,dstr] = parse_xy(dstr)
    [x,dstr] = parse_x(dstr);
    [y,dstr] = parse_x(dstr);
    if isempty(y) 
      xy = [];
      return;
    end
    xy = [x y];
  end
  function [key,dstr] = parse_key(dstr)
    % only support open cubic Bezier splines for now
    keys = 'CcHhLlMmSsVvZz';
    [key,count,~,pos] = sscanf(dstr,'%c',1);
    if count~= 1 
      key = [];
    elseif ~ismember(key,keys)
      key = [];
    else
      dstr = dstr(pos:end);
    end
  end

  [key,dstr] = parse_key(dstr);
  assert(key == 'M','First key must be M');
  [P,dstr] = parse_xy(dstr);
  mi = size(P,1);
  C = [];
  z_seen = false;
  % I guess L is some kind of default...
  prev_key = 'L';
  while ~isempty(dstr)
    [key,dstr] = parse_key(dstr);
    % A command letter may be eliminated if an identical command letter would
    % otherwise precede it; for instance, the following contains an unnecessary
    % second "L" command:
    if isempty(key) && ~isempty(prev_key)
      key = prev_key;
    end
    Pabs = P(end,:);
    switch key
    case {'C','c'}
      C = [C;size(P,1)+[0 1 2 3]];
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      if key == 'c'
        P(end-2:end,:) = P(end-2:end,:) + Pabs;
      end
    case {'L','l','V','v','H','h'}
      C = [C;size(P,1)+[0 1 2 3]];
      switch key
      case {'L','l'}
        [XY,dstr] = parse_xy(dstr);
        Pnext = (key=='l')*Pabs + XY;
      case {'V','v'}
        [Y,dstr] = parse_x(dstr);
        Pnext = [1 (key=='v')].*Pabs + [0 Y];
      case {'H','h'}
        [X,dstr] = parse_x(dstr);
        Pnext = [(key=='h') 1].*Pabs + [X 0];
      end
      P(end+1,:) = Pabs + (1/3)*(Pnext-Pabs);
      P(end+1,:) = Pabs + (2/3)*(Pnext-Pabs);
      P(end+1,:) = Pnext;
    case {'S','s'}
      C = [C;size(P,1)+[0 1 2 3]];
      P(end+1,:) = P(end,:)+ P(end,:)-P(end-1,:);
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      if key == 's'
        P(end-1:end,:) = P(end-1:end,:) + Pabs;
      end
    case {'M','m'}
      [Pnext,dstr] = parse_xy(dstr);
      P = [P;Pnext+(key=='m')*Pabs];
      mi = size(P,1);
    case {'Z','z'}
      %assert(isempty(dstr),'Z Should be last key');
      % gobble white space
      dstr = strip(dstr);
      %if all(P(end,:)==P(mi,:))
      % augh I hate that I'm using epsilon here. probably the interp1 above is
      % leading to small numerical noise.
      if sum((P(end,:)-P(mi,:)).^2)<eps
        % close up naturally by identifying first and last point
        C(end,4) = mi;
        P = P(1:end-1,:);
      end
    otherwise
      error('%c key not supported',key)
    end
    %key
    %clf;hold on;arrayfun(@(c) set(plot_cubic(P(C(c,:),:)),'Color','b'),1:size(C,1));hold off;
    %pause
    prev_key = key;
  end
    

end
