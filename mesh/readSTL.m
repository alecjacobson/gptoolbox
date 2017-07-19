function [V,F,N] = readSTL(filename,varargin)
  % READSTL read a triangle mesh from an .stl file.
  %
  % [V,F] = readSTL(filename)
  % [V,F,N] = readSTL(filename,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   filename  path to .stl file
  %   Optional:
  %     'JoinCorners' followed by whether to join perfectly match corners to
  %     form a combinatorially connected mesh {false}.
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of faces
  %   N  #F by 3 list of face normals
  % 

  % default values
  join_corners = false;
  variable_name2 = [1,2,3];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'JoinCorners'}, {'join_corners'});
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


  is_ascii = false;
  fid = fopen(filename, 'r');
  header =fread(fid,80,'uchar=>schar'); % Read file title
  header = char(header(:)');
  is_ascii = startsWith(lower(header),'solid');
  fclose(fid);

  if is_ascii
    fid = fopen(filename, 'r');
    % discard header line
    fgets(fid);
    % The prefixing space is important here.
    D = fscanf(fid,' facet normal %f %f %f outer loop vertex %f %f %f vertex %f %f %f vertex %f %f %f endloop endfacet ');
    D = reshape(D,12,[])';
    N = D(:,1:3);
    V = reshape(D(:,4:12)',3,[])';
    F = reshape(1:size(V,1),3,size(V,1)/3)';
    fclose(fid);
  else
    [FVX,FVY,FVZ] = stlread(filename);
    V = [FVX(:) FVY(:) FVZ(:)];
    F = reshape(1:size(V,1),3,size(V,1)/3)';
    N = [];
  end

  if join_corners
    [V,~,J] = remove_duplicate_vertices(V,0);
    % remap faces
    F = J(F);
  end

end
