function [P,C,I,F,S,W,D] = readSVG_cubics(filename)
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
      f = empty_val;
      return;
    end
    if isempty(match{1}) || strcmp(match{1}{1},'none') || match{1}{1}(1) ~= '#'
      f = nan(1,3);
      return;
    end
    assert(numel(match{1}{1})==7);
    f = hex2dec(reshape(match{1}{1}(2:7),2,3)')';
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

  function [Pi,Ci] = ellipse_to_path(cx,cy,rx,ry)
    % special fraction
    s =  1193/2160;
    Pi = [1,0;1,s;s,1;0,1;-s,1;-1,s;-1,0;-1,-s;-s,-1;0,-1;s,-1;1,-s;1,0];
    Ci = [1,2,3,4;4,5,6,7;7,8,9,10;10,11,12,13];
    Pi = Pi.*[rx ry]+[cx cy];
  end

  function [Pi,Ci] = parse_ellipse_as_spline(ellipse)
    cx = sscanf(char(ellipse.getAttribute('cx')),'%g');
    cy = sscanf(char(ellipse.getAttribute('cy')),'%g');
    rx = sscanf(char(ellipse.getAttribute('rx')),'%g');
    ry = sscanf(char(ellipse.getAttribute('ry')),'%g');
    [Pi,Ci] = ellipse_to_path(cx,cy,rx,ry);
  end
  function [Pi,Ci] = parse_circle_as_spline(circle)
    cx = sscanf(char(circle.getAttribute('cx')),'%g');
    cy = sscanf(char(circle.getAttribute('cy')),'%g');
    r = sscanf(char(circle.getAttribute('r')),'%g');
    if isempty(r)
      r = 0;
    end
    [Pi,Ci] = ellipse_to_path(cx,cy,r,r);
  end

  function Pi = parse_rect_as_polygon(rect)
    x = sscanf(char(rect.getAttribute('x')),'%g');
    y = sscanf(char(rect.getAttribute('y')),'%g');
    if isempty(x); x = 0; end
    if isempty(y); y = 0; end
    w = sscanf(char(rect.getAttribute('width')),'%g');
    h = sscanf(char(rect.getAttribute('height')),'%g');
    Pi = [x y;x+w y;x+w y+h;x y+h;x y];
  end

  function T = get_transform(kid)
    [T,foundT] = sscanf(char(kid.getAttribute('transform')),'matrix(%g %g %g %g %g %g)');
    if foundT == 0
      T = eye(2,3);
    else
      assert(foundT == 6)
      T = reshape(T,2,3);
    end
  end

  function Pi = parse_polygon(poly)
    Pi = reshape(sscanf(char(poly.getAttribute('points')),'%g,%g '),2,[])';
    % if end isn't exactly begin, make it so
    if any(Pi(1,:) ~= Pi(end,:))
      Pi(end+1,:) = Pi(1,:);
    end
  end

  function Pi = parse_polyline(pline)
    Pi = reshape(sscanf(char(pline.getAttribute('points')),'%g,%g '),2,[])';
  end

  function Pi = parse_line(line)
    get_scalar = @(line,key) sscanf(char(line.getAttribute(key)),'%g');
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
        % convert to cubic
        closed = all(Pi(1,:)==Pi(end,:));
        Pii = interp1(1:size(Pi,1),Pi,linspace(1,size(Pi,1),(size(Pi,1)-1)*3+1));
        Cii = [1 2 3 4]+3*(0:size(Pi,1)-2)';
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
      case 'g'
        % recurse...
        g_kids = kid.getChildNodes();
        [g_P,g_C,g_I,g_F,g_S,g_W,g_D,g_k] = process_kids(g_kids);
        C = [C;size(P,1)+g_C];
        P = [P;g_P];
        I = [I;k+g_I];
        F = [F;g_F];
        S = [S;g_S];
        W = [W;g_W];
        D = [D;g_D];
        k = k+g_k;
        continue;
      case '#text'
        % knowingly skip
        continue;
      otherwise
        fprintf('Skipping %s...\n',name);
        continue;
      end
      Si = get_color(kid,'stroke',nan(1,3));
      Fi = get_color(kid,'fill',[0 0 0]);
      Wi = get_scalar(kid,'stroke-miterlimit');
      Di = get_not_none(kid,'display');
      Ti = get_transform(kid);
      Pii = [Pii ones(size(Pii,1),1)]*Ti';

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
  kids = svg.getChildNodes();
  [P,C,I,F,S,W,D,k] = process_kids(kids);

end
