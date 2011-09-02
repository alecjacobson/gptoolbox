function writePOLY(varargin)
  % WRITEPOLY prints a vertices to a .poly file, with E connecting those vertices
  %
  % writePOLY(poly_file_name,poly_struct)
  % writePOLY(poly_file_name,V,E,H)
  %
  % Input
  %   poly_file_name:  name of output file as string (caution! will clobber
  %                    existing)
  %   poly:      struct array where each element contains fields:
  %                    x,y,hole
  %     OR
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge E
  %   H  #H by 2 list of hole positions
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: png2objandtga, png2poly
  %


  if(nargin == 4)
    poly_file_name = varargin{1};
    V = varargin{2};
    E = varargin{3};
    H = varargin{4};
  elseif(nargin == 2)
    poly_file_name = varargin{1};
    poly = varargin{2};
    [V,E,H] = poly2VEH(poly);
  else
    error('Wrong number of inputs');
  end
  % open file for writing
  poly_file_handle = fopen(poly_file_name,'w');

  % vertices section
  fprintf(poly_file_handle,'# vertices\n');
  fprintf(poly_file_handle,'%d 2 0 0\n', size(V,1));
  for j=1:size(V,1),
    fprintf(poly_file_handle,'%d %.17f %.17f\n',j,V(j,1),V(j,2));
  end

  % E section
  fprintf(poly_file_handle,'# E\n');
  fprintf(poly_file_handle,'%d 0\n', size(E,1));
  for j=1:size(E,1),
    fprintf(poly_file_handle,'%d %d %d\n',j,E(j,1),E(j,2));
  end

  % holes section
  if(exist('H','var'))
    fprintf(poly_file_handle,'# holes\n');
    fprintf(poly_file_handle,'%d 0\n', size(H,1));
    for j=1:size(H,1),
      fprintf(poly_file_handle,'%d %.17f %.17f\n',j,H(j,1),H(j,2));
    end
  else
    % mandatory number of holes lines
    fprintf(poly_file_handle,'# holes\n');
    fprintf(poly_file_handle,'0\n');
  end

  fprintf(poly_file_handle,'\n');
  fclose(poly_file_handle);

end
