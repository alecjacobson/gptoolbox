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
  %   C  #C by 2 list of control points of all polys
  %   E  #E by 2 list of edge indices of all polys into C
  %

  tic;
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
    Pi = [x y;x+w y;x+w y+h;x y+h;x y];
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


  C = [];
  E = [];
  for pi = 1:numel(P)
    E = [E;size(C,1)+[1:size(P{pi},1)-1;2:size(P{pi},1)]'];
    C = [C;P{pi}];
  end

end
