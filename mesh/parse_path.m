function [P,C,Pabs] = parse_path(dstr,varargin)
  % [P,C] = parse_path(dstr)
  %
  % Inputs:
  %   dstr  string of the value of the 'd' attribute
  % Outputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %   Pabs  2d position of cursor
  %
  % Known issues:
  %   This implementation is surely O(#dstr² + #P² + #C²). To fix this, we need
  %   to avoid using sscanf which appears to require copying to read from the
  %   middle of a string and switch P=[P;Pi] style concatenations to amortized
  %   P=[P Pi]; concatenations.
  %


  % These are all annoyingly O(#dstr) 
  function [x,dstr] = parse_x(dstr)
    [x,count,~,pos] = sscanf(dstr,'%g',1);
    if count~= 1
      x = [];
      return;
    end
    dstr = dstr(pos:end);
    dstr = strip(strip(dstr,'left',','),'left',' ');
  end
  function [x,dstr] = parse_flag(dstr)
    pos = 1;
    % eat delimiters
    while pos <= length(dstr) && ...
        (dstr(pos) == ' ' || dstr(pos) == ',')
      pos = pos + 1;
    end
    if dstr(pos) == '0' || dstr(pos) == '1'
      x = dstr(pos) == '1';
      pos = pos + 1;
    else
      x = [];
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
    keys = 'CcHhLlMmSsVvZzQqAaTt';
    [key,count,~,pos] = sscanf(dstr,'%c',1);
    if count~= 1 
      key = [];
    elseif ~ismember(key,keys)
      key = [];
    else
      dstr = dstr(pos:end);
    end
  end

  data = [];
  Pabs = [0 0];
  split = true;
  debug = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Cursor','Debug','Split'}, ...
    {'Pabs','debug','split'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end


  if split
    % This is band-aid for the problem that parse_path is O(n²). Really we
    % should be doing this for all commands not just Mm.
    %
    % Split into sub-paths
    pattern = '([Mm][^Mm]*)';
    S = regexp(dstr, pattern, 'match');
    P = [];
    C = [];
    Pabs = [0 0];
    for i = 1:length(S)
      [Pi,Ci,Pabs] = parse_path(S{i},'Cursor',Pabs,'Split',false);
      C = [C size(P,2)+Ci'];
      P = [P Pi'];
    end
    P = P';
    C = C';
    return;
  end


  dstr = strrep(dstr,',',' ');
  dstr = regexprep(dstr,'(&#x[A-z0-9];)+',' ');
  dstr = regexprep(dstr,'\s+',' ');

  if isempty(dstr)
    P = zeros(0,2);
    C = zeros(0,4);
    return;
  end

  [key,dstr] = parse_key(strtrim(dstr));
  assert(key == 'M' || key == 'm','First key must be M or m');
  [P,dstr] = parse_xy(dstr);
  P = (key == 'm')*Pabs + P;
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

      %% Handle full circles made of two 180° arcs as a special case.
      %if key == 'a'
      %end

      %[Z,count,~,pos] = sscanf(dstr,'%g',7);
      %assert(count == 7,'Expected 7 values for A,a command');
      %% O(n) :-(
      %dstr = strtrim(dstr(pos:end));
      [rx,dstr] = parse_x(dstr);
      [ry,dstr] = parse_x(dstr);
      [phi,dstr] = parse_x(dstr);
      % it seems valid for the flags to appear immediately adjacent to each
      % other like '10' or '11' which should be read as '1 0' and '1 1' rather
      % than 10→1 and 11→1
      % 
      % Another file ('222268.svg') appeared to use '1,11' to mean '1 1' but I
      % can't think of a way to correctly handle both uses.
      [large_arc,dstr] = parse_flag(dstr);
      [sweep,dstr] = parse_flag(dstr);
      [Pnext,dstr] = parse_xy(dstr);
      if isempty(Pnext)
        error('parse_path: wrong number of args for A,a');
      end

      Pnext = (key=='a')*Pabs + Pnext;

      if ~isequal(Pabs,Pnext)
        [Pe,Ce] = arc_to_cubics(Pabs,Pnext,rx,ry,phi,large_arc,sweep);

        if ~isempty(Ce)
          C = [C;size(P,1)+Ce-1];
          P = [P;Pe(2:end,:)];
          Pabs = P(end,:);
        end
      end
    case {'T','t'}
      Q1 = Pabs;
      if ismember(prev_key,'QqTt')
        Q2 = Pabs + (Pabs-Qprev);
      else
        Q2 = Pabs;
      end
      [Q3,dstr] = parse_xy(dstr);
      Q3 = Q3 + (key=='t')*Pabs;
      Q = [Q1;Q2;Q3];
      CQ = quadratic_to_cubic(Q);
      C = [C;size(P,1)+[0 1 2 3]];
      P = [P;CQ(2:end,:)];
      Qprev = Q(2,:);
      Pabs = P(end,:);

    case {'Q','q'}
      Q1 = Pabs;
      [Q2,dstr] = parse_xy(dstr);
      [Q3,dstr] = parse_xy(dstr);
      Q = [Q1;[Q2;Q3] + (key=='q')*Pabs];
      CQ = quadratic_to_cubic(Q);
      C = [C;size(P,1)+[0 1 2 3]];
      P = [P;CQ(2:end,:)];

      Qprev = Q(2,:);
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
