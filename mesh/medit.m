function [s,r] = medit(varargin)
  % MEDIT wrapper to launch medit to visualize a given tetmesh
  %
  % [s,r] = medit(V,T,F)
  % [s,r] = medit(V,T,F,'ParameterName',ParameterValue)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   T  #T by 4|5 list of tetrahedron findices (additional column is color
  %     index)
  %   F  #F by 3|4 list of face indices (additional column is color index)
  %   Optional:
  %     'Wait'  followed by boolean whether to wait for completion {true}
  %     'Data'  followed by #V|#T+#F data values {[]}
  % Outputs:
  %   s,r  result of system call
  %

  V = varargin{1};
  T = varargin{2};
  F = varargin{3};
  wait = true;
  D = [];

  ii = 4;
  while(ii<=nargin)
    switch varargin{ii}
    case 'Wait'
      ii = ii + 1;
      assert(ii<=nargin);
      wait = varargin{ii};
    case 'Data'
      ii = ii + 1;
      assert(ii<=nargin);
      D = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii+1;
  end

  if ~exist('wait','var')
    wait = true;
  end

  % Change these paths accordingly
  MEDIT_PATH = '/usr/local/igl/igl_lib/external/medit/medit';
  %MEDIT_PATH = '/opt/local/bin/medit';
  TEMP_MESH_FILE  = '/var/tmp/temp.mesh';
  TEMP_MEDIT_FILE = '/var/tmp/temp.medit';
  TEMP_BB_FILE = '/var/tmp/temp.bb';

  % write temporary mesh
  writeMESH(TEMP_MESH_FILE,V,T,F);

  % write default medit options
  f = fopen(TEMP_MEDIT_FILE,'w');
  fprintf(f, [...
    'BackgroundColor 1 1 1\n' ...
    'LineColor 0 0 0\n'       ...
    'WindowSize 1024 800\n']);
  if size(F,2) > 3 || size(T,2) > 4
    fprintf(f,['NbMaterials\n']);
    %% collect indices
    %I = [F(:,3:end);T(:,4:end)];
    %count = max(I) - min(I);
    %fprintf(f,'%d\n',count);
    fprintf(f,'RenderMode colorshadinglines\n');
  else
    fprintf(f,'RenderMode shading + lines\n');
  end
  fclose(f);

  % write bb file
  if isempty(D)
    if exist(TEMP_BB_FILE,'file')
      delete(TEMP_BB_FILE);
    end
  else
    switch size(D,1)
    case size(T,1)
      warning('Appending junk for face data since medit needs it.');
      D = [D;linspace(min(D),max(D),size(F,1))'];
      type = 1;
    case size(V,1)
      type = 2;
    case size(T,1)+size(F,1)
      type = 1;
    otherwise
      error('size(Data,1) %d should be == size(V,1) (%d) or size(T,1) + size(F,1) (%d + %d = %d)',size(D,1),size(V,1),size(T,1),size(F,1),size(T,1)+size(F,1));
    end
    writeBB(TEMP_BB_FILE,D,type);
  end

  command = [MEDIT_PATH ' ' TEMP_MESH_FILE ' ' TEMP_MEDIT_FILE];
  if ~wait
    command = [command ' &'];
  end
  [s,r] = system(command);
end
