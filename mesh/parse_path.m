function [P,C,state] = parse_path(dstr,varargin)
  % [P,C] = parse_path(dstr)
  %
  % Inputs:
  %   dstr  string of the value of the 'd' attribute
  % Outputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %   state
  %
  % Known issues:
  %   This implementation is surely O(#dstr² + #P² + #C²). To fix this, we need
  %   to avoid using sscanf which appears to require copying to read from the
  %   middle of a string and switch P=[P;Pi] style concatenations to amortized
  %   P=[P Pi]; concatenations.
  %


  % These should now be output sensitive (i.e., O(m) not O(n) where m is
  % position of final character needed to parse and n is length of dstr)
  function [x,pos] = parse_x(dstr,pos)
    [x,count,~,pos_off] = sscanf(dstr(pos:end),'%g',1);
    if count~= 1
      x = [];
      return;
    end
    pos = pos + pos_off - 1;
    % Eat spaces and commas
    while pos <= length(dstr) && (dstr(pos) == ',' || dstr(pos) == ' '); pos = pos + 1; end
  end
  function [x,pos] = parse_flag(dstr,pos)
    % Eat spaces and commas
    while pos <= length(dstr) && (dstr(pos) == ',' || dstr(pos) == ' '); pos = pos + 1; end
    if dstr(pos) == '0' || dstr(pos) == '1'
      x = dstr(pos) == '1';
      pos = pos + 1;
    else
      x = [];
    end
    % Eat spaces and commas
    while pos <= length(dstr) && (dstr(pos) == ',' || dstr(pos) == ' '); pos = pos + 1; end
  end
  function [xy,pos] = parse_xy(dstr,pos)
    [x,pos] = parse_x(dstr,pos);
    [y,pos] = parse_x(dstr,pos);
    if isempty(y) 
      xy = [];
      return;
    end
    xy = [x;y];
  end

  % Annoying way to effectively impement a class.
  function [state,pos] = parse_A(state,dstr,pos)
    [rx, pos] = parse_x(dstr,pos);
    [ry, pos] = parse_x(dstr,pos);
    [phi,pos] = parse_x(dstr,pos);
    % it seems valid for the flags to appear immediately adjacent to each
    % other like '10' or '11' which should be read as '1 0' and '1 1' rather
    % than 10→1 and 11→1
    % 
    % Another file ('222268.svg') appeared to use '1,11' to mean '1 1' but I
    % can't think of a way to correctly handle both uses.
    [large_arc,pos] = parse_flag(dstr,pos);
    [sweep,pos] = parse_flag(dstr,pos);
    [Pnext,pos] = parse_xy(dstr,pos);
    if isempty(Pnext)
      error('parse_path: wrong number of args for A,a');
    end

    Pnext = (state.key=='a')*state.Pabs + Pnext;

    if ~isequal(state.Pabs,Pnext)
      [Pe,Ce] = arc_to_cubics(state.Pabs',Pnext',rx,ry,phi,large_arc,sweep);

      if ~isempty(Ce)
        state.C = [state.C size(state.P,2)+Ce'-1];
        state.P = [state.P Pe(2:end,:)'];
        state.Pabs = state.P(:,end);
      end
    end
  end


  function [state,pos] = parse_C(state,dstr,pos)
    state.C = [state.C size(state.P,2)+[0 1 2 3]'];
    [xy,pos] = parse_xy(dstr,pos);
    state.P(:,end+1) = xy;
    [xy,pos] = parse_xy(dstr,pos);
    state.P(:,end+1) = xy;
    [xy,pos] = parse_xy(dstr,pos);
    state.P(:,end+1) = xy;

    if state.key == 'c'
      state.P(:,end-2:end) = state.P(:,end-2:end) + state.Pabs;
    end
    state.Pabs = state.P(:,end);
  end

  function [state,pos] = parse_L(state,dstr,pos)
      % augh I hate that I'm using epsilon here. probably the interp1 above is
      % leading to small numerical noise.
      if (state.key == 'Z' || state.key == 'z') && size(state.P,2)==state.mi
        % degenerate single point, ignore
      elseif (state.key == 'Z' || state.key == 'z') && sum((state.P(:,end)-state.P(:,state.mi)).^2)<eps
        % close up naturally by identifying first and last point
        if ~isempty(state.C)
          state.C(4,end) = state.mi;
        end
        % Pop last point
        state.P =    state.P(:,1:end-1);
        state.Pabs = state.P(:,state.mi);
      else
        state.C = [state.C size(state.P,2)+[0 1 2 3]'];
        switch state.key
        case {'L','l'}
          [XY,pos] = parse_xy(dstr,pos);
          Pnext = (state.key=='l')*state.Pabs + XY;
        case {'V','v'}
          [Y,pos] = parse_x(dstr,pos);
          Pnext = [1;(state.key=='v')].*state.Pabs + [0;Y];
        case {'H','h'}
          [X,pos] = parse_x(dstr,pos);
          Pnext = [(state.key=='h');1].*state.Pabs + [X;0];
        case {'Z','z'}
          % I don't see how Z or z are different.
          % Don't parse xy, connect to beginning 
          Pnext = state.P(:,state.mi);
        end
        state.P(:,end+1) = state.Pabs + (1/3)*(Pnext-state.Pabs);
        state.P(:,end+1) = state.Pabs + (2/3)*(Pnext-state.Pabs);
        state.P(:,end+1) = Pnext;
        state.Pabs = state.P(:,end);
      end
  end

  function [state,pos] = parse_M(state,dstr,pos)
    [Pnext,pos] = parse_xy(dstr,pos);
    Pm = Pnext+(state.key=='m')*state.Pabs;
    state.P = [state.P Pm];
    state.mi = size(state.P,2);
    state.Pabs = state.P(:,end);
    % Augh this can even happen after 'mz', set key here so below state.prev_key is
    % correct.
    state.key = char('L' + (state.key-'M'));
  end

  function [state,pos] = parse_S(state,dstr,pos)
    state.C = [state.C size(state.P,2)+[0 1 2 3]'];
    if ismember(state.prev_key,'SsCc')
      state.P(:,end+1) = state.P(:,end)+ state.P(:,end)-state.P(:,end-1);
    else
      state.P(:,end+1) = state.P(:,end);
    end
    [state.P(:,end+1),pos] = parse_xy(dstr,pos);
    [state.P(:,end+1),pos] = parse_xy(dstr,pos);
    if state.key == 's'
      state.P(:,end-1:end) = state.P(:,end-1:end) + state.Pabs;
    end
    state.Pabs = state.P(:,end);
  end

  function [state,pos] = parse_Q(state,dstr,pos)
      Q1 = state.Pabs;
      [Q2,pos] = parse_xy(dstr,pos);
      [Q3,pos] = parse_xy(dstr,pos);
      Q = [Q1 [Q2 Q3] + (state.key=='q')*state.Pabs]';
      CQ = quadratic_to_cubic(Q);
      state.C = [state.C size(state.P,2)+[0 1 2 3]'];
      state.P = [state.P CQ(2:end,:)'];

      state.Qprev = Q(2,:)';
      state.Pabs = state.P(:,end);
  end
  function [state,pos] = parse_T(state,dstr,pos)
      Q1 = state.Pabs;
      if ismember(state.prev_key,'QqTt')
        Q2 = state.Pabs + (state.Pabs-state.Qprev);
      else
        Q2 = state.Pabs;
      end
      [Q3,pos] = parse_xy(dstr,pos);
      Q3 = Q3 + (state.key=='t')*state.Pabs;
      Q = [Q1 Q2 Q3]';
      CQ = quadratic_to_cubic(Q);
      state.C = [state.C size(state.P,2)+[0 1 2 3]'];
      state.P = [state.P CQ(2:end,:)'];
      state.Qprev = Q(2,:);
      state.Pabs = state.P(:,end);
  end

  data = [];
  state_Pabs = [0 0];
  split = true;
  debug = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Cursor','Debug','Split'}, ...
    {'state_Pabs','debug','split'});
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
  state = struct();
  state.Pabs = reshape(state_Pabs,[],1);


  if split
    % This is band-aid for the problem that parse_path is O(n²). Really we
    % should be doing this for all commands not just Mm.
    %
    % Split into sub-paths
    pattern = '([Mm][^Mm]*)';
    S = regexp(dstr, pattern, 'match');
    state.P = [];
    state.C = [];
    state.Pabs = [0 0];
    for i = 1:length(S)
      [Pi,Ci,state_i] = parse_path(S{i},'Cursor',state.Pabs,'Split',false);
      state.Pabs = state_i.Pabs;
      state.C = [state.C size(state.P,2)+Ci'];
      state.P = [state.P Pi'];
    end
    P = state.P';
    C = state.C';
    return;
  end

  % Combine these somehow?
  dstr = strrep(dstr,',',' ');
  dstr = regexprep(dstr,'(&#x[A-z0-9];)+',' ');
  dstr = regexprep(dstr,'\s+',' ');

  if isempty(dstr)
    P = zeros(0,2);
    C = zeros(0,4);
    return;
  end

  keys = 'CcHhLlMmSsVvZzQqAaTt';
  pattern = sprintf('([%s][^%s]*)',keys,keys);
  S = regexp(dstr, pattern, 'match');
  if isempty(S)
    P = zeros(0,2);
    C = zeros(0,4);
    return;
  end

  state.P = zeros(2,0);
  state.C = zeros(4,0);
  state.key = S{1}(1);
  state.Qprev = [];
  assert(state.key == 'M' || state.key == 'm','First key must be M or m');
  for si = 1:numel(S)
    dstr_i = S{si};
    state.key = dstr_i(1);
    pos = 2;
    % Since these are each copying and editing dstr_i, this still incurs a O(n²)
    % performance hit, where n is the length of a command (or chain of commands
    % with same key [potentially very long]).
    %
    % A finally optimization would be to be careful to not copy/edit dstr_i or
    % to do an initial pass to break up chains of commands.
    %
    % Currently, I think it would be easier to maintain dstr_i and pos and
    % change each parse_ command to increment on pos without altering dstr_i
    while true
      switch(state.key)
      case {'A','a'}
        [state,pos] = parse_A(state,dstr_i,pos);
      case {'C','c'}
        [state,pos] = parse_C(state,dstr_i,pos);
      case {'L','l','V','v','H','h','Z','z'}
        [state,pos] = parse_L(state,dstr_i,pos);
      case {'M','m'}
        [state,pos] = parse_M(state,dstr_i,pos);
      case {'Q','q'}
        [state,pos] = parse_Q(state,dstr_i,pos);
      case {'S','s'}
        [state,pos] = parse_S(state,dstr_i,pos);
      case {'T','t'}
        [state,pos] = parse_T(state,dstr_i,pos);
      otherwise
        error('%c key not supported',state.key)
      end

      % would like to get rid of this:
      while pos <= length(dstr_i) && (dstr_i(pos) == ',' || dstr_i(pos) == ' '); pos = pos + 1; end
      if pos>length(dstr_i)
        break;
      end
    end
    state.prev_key = state.key;
  end


  P = state.P';
  C = state.C';
end
