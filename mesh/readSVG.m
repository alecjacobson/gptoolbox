function [P,S,C,E] = readSVG(filename)
  % READSVG  Read <polygon>, <polyline> from a .svg file
  %
  % P = readSVG(filename)
  %
  % Input:
  %   filename  path to .svg file
  % Output:
  %   P  #P list of polylines, if the polygon is closed the the last point will
  %      be the same as the first.
  %   S  #P by 3 list of stroke colors
  %   C  #C list of cubic bezir paths
  %
  % Note: order is **not** maintained
  %
  % Example:
  %   % read polylines from an svg 
  %   [P,S] = readSVG('polylines.svg');

  tic;
  % xmlread appears to be buggy if filename contains ~ etc.
  if isunix
    [~,filename] = system(sprintf('find %s',filename));
    filename = strtrim(filename);
  end
  xDoc = xmlread(filename);
  toc
  
  P = {};
  S = [];

  rects = xDoc.getElementsByTagName('rect');
  for pi = 0:rects.getLength-1
    rect = rects.item(pi);
    x = sscanf(char(rect.getAttribute('x')),'%g');
    y = sscanf(char(rect.getAttribute('y')),'%g');
    w = sscanf(char(rect.getAttribute('width')),'%g');
    h = sscanf(char(rect.getAttribute('height')),'%g');
    [T,foundT] = sscanf(char(rect.getAttribute('transform')),'matrix(%g %g %g %g %g %g)');
    if foundT == 0
      T = eye(2,3);
    else
      assert(foundT == 6)
      T = reshape(T,2,3);
    end
    Pi = [[x y;x+w y;x+w y+h;x y+h;x y] ones(5,1)]*T';
    P{end+1} = Pi;
    hex = char(rect.getAttribute('stroke'));
    if ~isempty(hex)
      S(end+1,:) = hex2rgb(hex(2:end));
    end
  end

  % Read polygons
  polys = xDoc.getElementsByTagName('polygon');
  for pi = 0:polys.getLength-1
    poly = polys.item(pi);
    Pi = reshape(sscanf(char(poly.getAttribute('points')),'%g,%g '),2,[])';
    % if end isn't exactly begin, make it so
    if any(Pi(1,:) ~= Pi(end,:))
      Pi(end+1,:) = Pi(1,:);
    end
    P{end+1} = Pi;
    hex = char(poly.getAttribute('stroke'));
    if ~isempty(hex)
      S(end+1,:) = hex2rgb(hex(2:end));
    end
  end

  plines = xDoc.getElementsByTagName('polyline');
  for pi = 0:plines.getLength-1
    pline = plines.item(pi);
    P{end+1} = reshape(sscanf(char(pline.getAttribute('points')),'%g,%g '),2,[])';
    hex = char(pline.getAttribute('stroke'));
    if ~isempty(hex)
      S(end+1,:) = hex2rgb(hex(2:end));
    end
  end

  % Combine all into segment complex
  C = [];
  E = [];
  for ii = 1:numel(P)
    E = [E;size(C,1)+[1:size(P{ii},1)-1;2:size(P{ii},1)]'];
    C = [C;P{ii}];
  end
    
  


end
