function writeBF(filename,V,WI,P)
  % WRITEBF
  %
  % writeBF(filename,V,WI,P)
  %
  % Write a bone forest to a .bf file
  %
  % Input:
  %  filename  .bf file name
  %  V  # vertices by 3 list of vertex offsets from parents
  %  WI  # vertices list of indices of weights (0 means no weights)
  %  P  # vertices list of indices of bone parents (0 means root)
  % 

  if ~exist('WI','var')
    WI = 1:size(V,1);
    warning('Weight indices not given, assuming weight indices correspond to order');
  end

  if ~exist('P','var')
    P = zeros(size(V,1),1);
    warning('Parent indices not given, assuming all are roots');
  end

  if size(V,2) == 2
    warning('Vertex positions given in 2D, assuming 0 for all z-coordinates');
    V = [V zeros(size(V,1),1)];
  end

  assert(numel(P) == size(V,1));
  assert(numel(WI) == size(V,1));

  % make WI and P into columns to match V
  WI = WI(:);
  P = P(:);

  % P and WI are stored in matlab style one-indexing but should be printed in
  % c++ style zero-indexing
  P = P-1;
  WI = WI-1;

  disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  fprintf(fp,'%d %d %g %g %g\n',[WI';P';V']);
  fclose(fp);
end
