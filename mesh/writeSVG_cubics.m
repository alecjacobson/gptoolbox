function writeSVG_cubics(filename,P,C,I,F,S,W,units)
  if ~exist('units','var')
    units = 'mm';
  end
  fh = fopen(filename,'w');
  fprintf(fh,['<svg version="1.1" id="Layer_1" ' ...
    ' xmlns="http://www.w3.org/2000/svg" ' ...
    ' xmlns:xlink="http://www.w3.org/1999/xlink" ' ...
    ' width="%d%s" height="%d%s" ' ...
    ' viewBox="%d %d %d %d" ' ...
    ' xml:space="preserve">\n'], range(P(:,1)),units,range(P(:,2)),units,min(P),range(P));
  w = min(0.001*normrow(max(P)-min(P)),1);
  for c = 1:size(C,1)
    fprintf(fh,'<path ');
    fprintf(fh,'style="fill:none;stroke:#000000;stroke-width:%f;stroke-miterlimit:10;" ',w);
    fprintf(fh,'d="M%f,%fC%f,%f,%f,%f,%f,%f"></path>\n',P(C(c,:),:)');
  end
  fprintf(fh,'</svg>\n');
  fclose(fh);
end
