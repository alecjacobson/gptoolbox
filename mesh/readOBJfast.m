function [V,F] = readOBJfast(filename,varargin)
  % readOBJfast
  %
  % reads an OBJ file quickly, but OBJ file should be formated *simply*
  %
  % Input:
  %  filename  path to .obj file
  % Outputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %
  % filename should be a text file containing first vertex positions where each
  % line is:
  % v x y z
  % then a list of faces where each line is:
  % f i j k
  % and exactly that. Any lines before vertices, or after vertices but before
  % faces are ignored. Comments throughout vetex and face lines will surely
  % break this.
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % For complete(r) support use readOBJ
  % 
  % See also readOBJ
  %

  function file_wrangler()
    error_handler = @final_error;
    if isunix()
      % try to make the file easier to read
      tmpf = tempname();
      escaped_filename = strrep(filename,' ','\ ');
      cmd = sprintf( ...
        ['echo -e "$(grep "^v " %s) \\n$(grep "^f " %s | ' ...
        'sed -e "s/\\([0-9]*\\)\\(\\/[0-9]*\\)*/\\1/g")" >%s'], ...
        escaped_filename,escaped_filename,tmpf);
      s = system(cmd);
      if s==0
        readOBJfast_helper(tmpf);
        %delete(tmpf);
        return;
      end
    end
    final_error();
  end

  function final_error()
    error('Bad format. Try readOBJ...\n');
  end

  function readOBJfast_helper(filename)
    fp = fopen(filename);
    % read in vertices
    V = [];
    while(true)
      V = fscanf(fp,' v %lg %lg %lg ',inf);
      if(prod(size(V)) > 0)
        break;
      else
        line = fgets(fp);
        if(prod(size(line)) == 0)
          fclose(fp);
          error_handler();
          return;
        end
      end
    end
    V = reshape(V,3,size(V,1)/3)';
    % read in faces
    F = [];
    if quads
      format = ' f %d %d %d %d ';
    else
      format = ' f %d %d %d ';
    end
    while(true)
      F = fscanf(fp,format,inf);
      if(numel(F) > 0)
        break;
      else
        line = fgets(fp);
        if(numel(line) == 0 || (isnumeric(line) && line == -1))
          fclose(fp);
          error_handler();
          return;
        end
      end
    end
    if quads
      if mod(size(F,1),4) ~= 0
        fclose(fp);
        error_handler();
        return;
      end
      F = reshape(F,4,size(F,1)/4)';
    else
      if mod(size(F,1),3) ~= 0
        fclose(fp);
        error_handler();
        return;
      end
      F = reshape(F,3,size(F,1)/3)';
    end
    if ~feof(fp)
      fclose(fp);
      error_handler();
      return;
    end
    fclose(fp);
  end

  % default values
  quads = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Quads'}, ...
    {'quads'});
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

  error_handler = @file_wrangler;
  readOBJfast_helper(filename)
end
