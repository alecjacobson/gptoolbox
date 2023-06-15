function [P,C] = parse_path(dstr)
  % [P,C] = parse_path(dstr)
  %
  % Inputs:
  %   dstr  string of the value of the 'd' attribute
  % Outputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %

  function [c,r] = circle_center(A, B, r, sweep)
    % Compute the midpoint M of A and B
    M = (A + B) / 2;
    
    % Compute the distance between A and M
    d = norm(A - M);
    
    % Calculate the distance from M to the center C
    if (r^2 - d^2)<0
      warning('A,a radius too small');
      %r = d+eps(abs(d));
      r = d+1e-7;
    end
    h = sqrt(r^2 - d^2);
    
    % Determine the direction from A to B
    dir = (B - A) / norm(B - A);
    
    % Determine the direction perpendicular to AB
    perp_dir = [-dir(2), dir(1)];
    
    % Compute the two possible centers
    c = M + ((sweep==1)*1+(sweep~=1)*-1)  * h * perp_dir;

end

  function [x,dstr] = parse_x(dstr)
    [x,count,~,pos] = sscanf(dstr,'%g',1);
    if count~= 1
      x = [];
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
    keys = 'CcHhLlMmSsVvZzQqAa';
    [key,count,~,pos] = sscanf(dstr,'%c',1);
    if count~= 1 
      key = [];
    elseif ~ismember(key,keys)
      key = [];
    else
      dstr = dstr(pos:end);
    end
  end

  dstr = strrep(dstr,'\n',' ');
  dstr = strrep(dstr,',',' ');

  [key,dstr] = parse_key(strtrim(dstr));
  assert(key == 'M' || key == 'm','First key must be M or m');
  [P,dstr] = parse_xy(dstr);
  mi = size(P,1);
  % Pabs will track the current cursor position
  Pabs = P(end,:);
  C = [];
  z_seen = false;
  % I guess L is some kind of default...
  prev_key = char('L' + (key-'M'));
  while ~isempty(dstr)
    [key,dstr] = parse_key(dstr);
    % A command letter may be eliminated if an identical command letter would
    % otherwise precede it; for instance, the following contains an unnecessary
    % second "L" command:
    if isempty(key) && ~isempty(prev_key)
      key = prev_key;
    end
    switch key
    case {'A','a'}
      [Z,count,~,pos] = sscanf(dstr,'%g',7);
      assert(count == 7,'Expected 7 values for A,a command');
      % O(n) :-(
      dstr = strtrim(dstr(pos:end));
      rx = Z(1);
      ry = Z(2);
      % x-axis rotation
      phi = Z(3);
      % large-arc-flag
      large_arc = Z(4);
      % sweep-flag
      sweep = Z(5);
      assert(phi == 0,'phi ~= 0 not supported');
      Pnext = (key=='a')*Pabs + [Z(6) Z(7)];

      assert(rx == ry,sprintf('rx (%g) ~= ry (%g) not supported',rx,ry));
      r = rx;
      [c,r] = circle_center(Pabs, Pnext, r, xor(sweep,large_arc));

      [Pe,Ce] = ellipse_to_spline(c(1),c(2),r,r);

      is_flat_tol = 0.01;
      Pcut = [Pabs;Pnext];
      J = zeros(size(Pe,1),1);
      for p = 1:size(Pcut,1)
        [~,Ie,Se,K] = point_spline_squared_distance(Pcut(p,:),Pe,Ce,is_flat_tol);
        if Se == 0
          J(Ce(Ie,1)) = p;
        elseif Se == 1
          J(Ce(Ie,4)) = p;
        else
          [Q1,Q2] = cubic_split(Pe(Ce(Ie,:),:),Se);
          J = [J(1:Ce(Ie,1));0;0;p;0;0;J(Ce(Ie,4):end)];
          Pe = [Pe(1:Ce(Ie,1),:);Q1(2:end,:);Q2(2:end-1,:);;Pe(Ce(Ie,4):end,:)];
          Ce = mod(((1:3:size(Pe,1)-1)'  + (0:3)) - 1,size(Pe,1))+1;
        end
      end
      I1 = find(J==1);
      I2 = find(J==2);

      % explicitly close up
      if isequal(Pe(1,:),Pe(end,:))
        Pe = Pe(1:end-1,:);
        J(1) = max(J(1),J(end));
      end
      Ce = mod(((1:3:size(Pe,1)-1)'  + (0:3)) - 1,size(Pe,1))+1;

      l = edge_lengths(Pe,Ce(:,[1 4]));

      % Mid points
      BC = ...
        (0.5).^3.*Pe(Ce(:,1),:) + ...
        3.*(0.5).^3.*Pe(Ce(:,2),:) + ...
        3.*(0.5).^3.*Pe(Ce(:,3),:) + ...
        (0.5).^3.*Pe(Ce(:,4),:);
      larger = sum((BC-Pabs).*((Pnext-Pabs)*[0 1;-1 0]),2)>0;
      if sum(l(larger)) < sum(l(~larger))
        larger = ~larger;
      end

      Ce = Ce(larger == large_arc,:);

      if ~any(Ce(:,1)==I1)
        Ce = fliplr(Ce);
      end
      assert( sum(Ce(:,1)==I1) == 1 )
      E = [Ce(:,1:2);Ce(:,2:3);Ce(:,3:4)];
      A = sparse(E(:,1),E(:,2),1,size(Pe,1),size(Pe,1));
      I = [I1];
      while I(end) ~= I2
        I = [I;find(A(I(end),:))];
      end
      Pe = Pe(I,:);
      Ce = (1:3:size(Pe,1)-1)'  + (0:3);
      C = [C;size(P,1)+Ce-1];
      P = [P;Pe(2:end,:)];

      Pabs = P(end,:);
    case {'C','c'}
      C = [C;size(P,1)+[0 1 2 3]];
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      if key == 'c'
        P(end-2:end,:) = P(end-2:end,:) + Pabs;
      end
      Pabs = P(end,:);
    case {'L','l','V','v','H','h','Z','z'}
      if (key == 'Z' || key == 'z')
        %assert(isempty(dstr),'Z Should be last key');
        % gobble white space
        dstr = strip(dstr);
      end
      % augh I hate that I'm using epsilon here. probably the interp1 above is
      % leading to small numerical noise.
      if (key == 'Z' || key == 'z') && size(P,1)==mi
        % degenerate single point, ignore
      elseif (key == 'Z' || key == 'z') && sum((P(end,:)-P(mi,:)).^2)<eps
        % close up naturally by identifying first and last point
        if ~isempty(C)
          C(end,4) = mi;
        end
        % Pop last point
        P = P(1:end-1,:);
        Pabs = P(mi,:);
      else
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
        case {'Z','z'}
          % I don't see how Z or z are different.
          % Don't parse xy, connect to beginning 
          Pnext = P(mi,:);
        end
        P(end+1,:) = Pabs + (1/3)*(Pnext-Pabs);
        P(end+1,:) = Pabs + (2/3)*(Pnext-Pabs);
        P(end+1,:) = Pnext;
        Pabs = P(end,:);
      end
    case {'S','s'}
      C = [C;size(P,1)+[0 1 2 3]];
      if ismember(prev_key,'SsCc')
        P(end+1,:) = P(end,:)+ P(end,:)-P(end-1,:);
      else
        P(end+1,:) = P(end,:);
      end
      [P(end+1,:),dstr] = parse_xy(dstr);
      [P(end+1,:),dstr] = parse_xy(dstr);
      if key == 's'
        P(end-1:end,:) = P(end-1:end,:) + Pabs;
      end
      Pabs = P(end,:);
    case {'M','m'}
      [Pnext,dstr] = parse_xy(dstr);
      Pm = Pnext+(key=='m')*Pabs;
      P = [P;Pm];
      mi = size(P,1);
      Pabs = P(end,:);
      % Augh this can even happen after 'mz', set key here so below prev_key is
      % correct.
      key = char('L' + (key-'M'));
    otherwise
      error('%c key not supported',key)
    end
    %key
    %clf;hold on;arrayfun(@(c) set(plot_cubic(P(C(c,:),:)),'Color','b'),1:size(C,1));hold off;
    %pause
    prev_key = key;
  end
    

end
