function [P,C,I,F,S,W,D] = readSVG_cubics(filename,varargin)
  % READSVG_CUBICS Read in paths and shapes from a .svg file, but convert
  % everything to cubic Bezier curves.
  %
  % [P,C,I,F,S,W,D] = readSVG_cubics(filename)
  % 
  % Inputs:
  %   filename  path to .svg file
  % Outputs:
  %   P  #P by 2 list of control point locations
  %   C  #C by 4 list of indices into P of cubics
  %   I  #C list of indices into total number of svg objects
  %   F  #max(I) list of fill colors, NaN means none
  %   S  #max(I) list of stroke colors, NaN means none
  %   W  #max(I) list of widths colors, NaN means none
  %   D  #max(I) list of display, 0 means 'none', 1 otherwise
  %
  % Example:
  %   clf;
  %   hold on;
  %   for c = 1:size(C,1)
  %   plot_cubic(P(C(c,:),:));
  %   end
  %   hold off
  %   set(gca,'YDir','reverse')
  %
  % Example:
  %   h=5;
  %   [P,C] = readSVG_cubics('~/Downloads/r.svg');
  %   [V,E] = spline_to_poly(P,C,0.01);
  %   V(:,2)=max(V(:,2))-V(:,2);
  %   [V,F] = triangle(V,E,[],'Flags',sprintf('-q33 -a%0.17f',median(edge_lengths(V,E))^2));
  %   [V,F] = extrude(V,F,'Levels',ceil(h/mean(mean(edge_lengths(V,F)))));
  %   V(:,3)=V(:,3)*h;
  %   tsurf(F,V);
  %   axis equal;
  %   writeOBJ('~/Repos/shape-editing-examples/r.obj',V,F);
  %   
  % Known issues:
  %    <use>, etc. not supported see below in code
  %    likely O(n²) in number of kids and O(m²) in the number of cubics
  %      (parse_path)


  function f = get_not_none(kid,key)
    style = char(kid.getAttribute('style'));
    [~,~,~,~,match] = regexp(style,[key ': *([^;]*);']);
    if ~isempty(match) && ~isempty(match{1}) && strcmp(match{1}{1},'none')
      f = false;
      return;
    end
    f = true;
  end

  function f = get_color(kid,key,empty_val)
    style = char(kid.getAttribute('style'));
    [~,~,~,~,match] = regexp(style,[key ': *([^;]*);']);
    if isempty(match)
      match = char(kid.getAttribute(key));
      if isempty(match)
        f = empty_val;
        return;
      end
      match = {{match}};
    end
    if isempty(match{1}) || strcmp(match{1}{1},'none') || match{1}{1}(1) ~= '#'
      f = nan(1,3);
      return;
    end
    assert(numel(match{1}{1})==4 || numel(match{1}{1})==7);
    f = hex2rgb(match{1}{1}(2:end));
    assert(all(size(f)==[1 3]));
  end
  function s = get_scalar(kid,key)
    [~,~,~,~,match] = regexp( ...
      char(kid.getAttribute('style')),[key ': *([^;]*);']);
    if isempty(match) || isempty(match{1}) || strcmp(match{1}{1},'none') 
      s = nan(1,1);
      return;
    end
    assert(numel(match{1}{1})>=1);
    s = str2num(match{1}{1});
  end


  function [Pi,Ci] = parse_ellipse_as_spline(ellipse)
    cx = sscanf(char(ellipse.getAttribute('cx')),'%g');
    cy = sscanf(char(ellipse.getAttribute('cy')),'%g');
    rx = sscanf(char(ellipse.getAttribute('rx')),'%g');
    ry = sscanf(char(ellipse.getAttribute('ry')),'%g');
    if isempty(cx)
      cx = 0; 
    end
    if isempty(cy)
      cy = 0; 
    end
    if isempty(rx)
      rx = 0; 
    end
    if isempty(ry)
      ry = 0; 
    end
    [Pi,Ci] = ellipse_to_spline(cx,cy,rx,ry);
  end
  function [Pi,Ci] = parse_circle_as_spline(circle)
    cx = sscanf(char(circle.getAttribute('cx')),'%g');
    cy = sscanf(char(circle.getAttribute('cy')),'%g');
    if isempty(cx)
      cx = 0; 
    end
    if isempty(cy)
      cy = 0; 
    end
    r = sscanf(char(circle.getAttribute('r')),'%g');
    if isempty(r)
      r = 0;
    end
    [Pi,Ci] = ellipse_to_spline(cx,cy,r,r);
  end

  function Pi = parse_rect_as_polygon(rect)
    x = sscanf(char(rect.getAttribute('x')),'%g');
    y = sscanf(char(rect.getAttribute('y')),'%g');
    if isempty(x); x = 0; end
    if isempty(y); y = 0; end
    w = sscanf(char(rect.getAttribute('width')),'%g');
    h = sscanf(char(rect.getAttribute('height')),'%g');
    if isempty(w) || isempty(h)
      warning('readSVG_cubics: auto rect.width or rect.height not supported');
      Pi = zeros(0,2);
      return;
    end
    if rect.hasAttribute('rx')
      warning('readSVG_cubics:rect.rx not supported');
    end
    if rect.hasAttribute('ry')
      warning('readSVG_cubics:rect.ry not supported');
    end
    Pi = [x y;x+w y;x+w y+h;x y+h;x y];
  end

  function T = get_transform(kid)
    T = parse_svg_transform(char(kid.getAttribute('transform')));
  end

  function Pi = parse_points_attribute(elem)
    points_str = char(elem.getAttribute('points'));
    points_str = strrep(points_str,',',' ');
    raw_numbers = sscanf(points_str,'%g');
    assert(mod(numel(raw_numbers),2)==0,'points attribute must have even number of numbers');
    Pi = reshape(raw_numbers,2,[])';
  end

  function Pi = parse_polygon(poly)
    Pi = parse_points_attribute(poly);
    % if end isn't exactly begin, make it so
    if any(Pi(1,:) ~= Pi(end,:))
      Pi(end+1,:) = Pi(1,:);
    end
  end

  function Pi = parse_polyline(pline)
    Pi = parse_points_attribute(pline);
  end

  function Pi = parse_line(line)
    function s = get_scalar(line,key)
      if line.hasAttribute(key)
        s = sscanf(char(line.getAttribute(key)),'%g');
      else
        % I guess 0 is default
        s = 0;
      end
    end
    Pi = [ ...
      get_scalar(line,'x1') get_scalar(line,'y1') ; ...
      get_scalar(line,'x2') get_scalar(line,'y2')];
  end

  function [P,C,I,F,S,W,D,k] = process_kids(kids)
    P = [];
    C = [];
    I = [];
    F = [];
    S = [];
    W = [];
    D = [];
    % number of actual svg objects encountered
    k = 0;
    for ii = 0:kids.getLength-1
      Pii = [];
      Cii = [];
      kid = kids.item(ii);
      name = char(kid.getNodeName());
      % strip 'svg:' if present as prefix
      if numel(name)>4 && strcmp(name(1:4),'svg:')
        name = name(5:end);
      end
      switch name
        % straight things
      case {'line','polyline','rect','polygon'}
        switch name
        case 'line'
          Pi = parse_line(kid);
        case 'polyline'
          Pi = parse_polyline(kid);
        case 'polygon'
          Pi = parse_polygon(kid);
        case 'rect'
          Pi = parse_rect_as_polygon(kid);
        otherwise
          fprintf('skipping %s...\n',name);
        end
        if size(Pi,1)<=1
          % skip
          Pii = zeros(0,2);
          Cii = zeros(0,4);
        else
          % convert to cubic
          Pii = interp1(1:size(Pi,1),Pi,linspace(1,size(Pi,1),(size(Pi,1)-1)*3+1));
          Cii = [1 2 3 4]+3*(0:size(Pi,1)-2)';
        end
        %if ~closed && filled
        %  % close it with a straight curve
        %  fprintf('should be closing...')
        %end
      case 'ellipse'
        [Pii,Cii] = parse_ellipse_as_spline(kid);
      case 'circle'
        [Pii,Cii] = parse_circle_as_spline(kid);
      case 'path'
        % I'm assuming path will not actually store something straight
        d = char(kid.getAttribute('d'));
        %if k == 6
        %  keyboard
        %end
        [Pii,Cii] = parse_path(d);
        %clf;hold on;arrayfun(@(c) set(plot_cubic(Pii(Cii(c,:),:)),'Color','b'),1:size(Cii,1));hold off;
        %pause
      case {'g','defs','clipPath'}
        % Illustrator appears to put clip-paths in <defs> at the top of the
        % file. Recurce into these similar to a group <g>
        g_kids = kid.getChildNodes();
        [g_P,g_C,g_I,g_F,g_S,g_W,g_D,g_k] = process_kids(g_kids);
        g_T = get_transform(kid);
        if ~isempty(g_T) && ~isempty(g_P)
          g_P = [g_P ones(size(g_P,1),1)]*g_T';
        end
        % augh. This is O(n²), should use a cell of transpose-transpose
        C = [C;size(P,1)+g_C];
        P = [P;g_P];
        I = [I;k+g_I];
        F = [F;g_F];
        S = [S;g_S];
        W = [W;g_W];
        D = [D;g_D];
        k = k+g_k;
        % I have a feeling that transforms applied to groups are not being read
        % correctly. See below. Seems like transforms are only applied to
        % Pii,Cii type of additions. So maybe groups should be writing into
        % Pii,Cii etc. and then added to P,C etc. below, rather than this
        % `continue`
        continue;
      case {'#comment','head','script','desc','text','#text','flowRoot','pattern','metadata','title','style','filter','linearGradient','radialGradient'}
        % intentionally skip
        continue;
      case {'use','font','marker','symbol','mask'}
        alert('readSVG_cubics:skip Skipping <%s> ...\n',name);
        % knowingly skip
        continue;
      otherwise
        % if name starts with sodipodi: or inkscape: 
        if startsWith(name,'sodipodi:') || startsWith(name,'inkscape:')
          % knowingly skip
          continue;
        end
        alert('readSVG_cubics:unknown Skipping unknown element <%s> ...\n',name);
        continue;
      end
      Si = get_color(kid,'stroke',nan(1,3));
      Fi = get_color(kid,'fill',[0 0 0]);

      Wi = get_scalar(kid,'stroke-miterlimit');
      Di = get_not_none(kid,'display');
      Ti = get_transform(kid);
      if ~isempty(Pii)
        Pii = [Pii ones(size(Pii,1),1)]*Ti';
      end

      % augh. This is O(n²), should use a cell of transpose-transpose
      if ~isempty(Cii)
        k = k+1;
        S = [S;Si];
        F = [F;Fi];
        W = [W;Wi];
        D = [D;Di];
        I = [I;repmat(k,size(Cii,1),1)];
        C = [C;size(P,1)+Cii];
        P = [P;Pii];
      end
    end
  end

  alert_type = 'error';

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'UnsupportedAlert'}, ...
    {'alert_type'});
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

  switch alert_type
  case 'error'
    alert = @(varargin) error(varargin{:}); 
  case 'warning'
    alert = @(varargin) warning(varargin{:}); 
  case 'ignore'
    alert = @(varargin) [];
  end


  % xmlread appears to be buggy if filename contains ~ etc.
  if isunix
    [~,filename] = system(sprintf('find %s',filename));
    filename = strtrim(filename);
  end
  xDoc = xmlread(filename);

  %paths = xDoc.getElementsByTagName('path');
  %P = [];
  %C = [];
  %for ii = 0:paths.getLength-1
  %  path = paths.item(ii);
  %  d = char(path.getAttribute('d'));
  %  [Pii,Cii] = parse_path(d);
  %  C = [C;size(P,1)+Cii];
  %  P = [P;Pii];
  %end

  svg = xDoc.getElementsByTagName('svg').item(0);
  if isempty(svg)
    % seems sometimes all of the tags are svg:svg, svg:path, … instead of svg,
    % path, … 
    svg = xDoc.getElementsByTagName('svg:svg').item(0);
  end
  kids = svg.getChildNodes();
  [P,C,I,F,S,W,D,k] = process_kids(kids);

end
