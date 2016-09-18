function C = depends(f,depth)
  % DEPENDS Computes dependencies for a given function, f.
  %
  % C = depends(f)
  % C = depends(f,depth)
  %
  % Inputs:
  %   f  string, name of function or path to function
  %   depth  number of levels to explore {Inf}
  % Outputs:
  %   C  cell array of paths to dependencies
  %
  % Example:
  %  % Get dependencies and zip into a file
  %  C = depends('myfun');
  %  zip('myfun.zip',C);
  %
  %  % Get dependencies and extract just file base names and print
  %  C = depends('myfun');
  %  N = regexprep(C,'^.*\/','')
  %  fprintf('myfun depends on:\n');
  %  fprintf('  %s\n',N{:})
  %
  %  % Get dependencies and remove mosek paths
  %  C = depends('myfun');
  %  C = C(cellfun(@isempty,strfind(C,'opt/local/mosek')));
  % 
  
  w = warning('off','MATLAB:DEPFUN:DeprecatedAPI');
  

  if(~exist('depth','var'))
    depth = Inf;
  end

  level = 0;
  C = {};
  C{end+1} = which(f);

  Q = {};
  Q = C;
  while( level < depth && ~isempty(Q))
    new_deps = {};
    while(~isempty(Q))
      p = Q{1};
      Q = Q(2:end);
      %new_deps = cat(1,new_deps,depfun(p,'-toponly','-quiet'));
      % This is the non-obsolete version but it's 100x slower : - (
      new_deps = cat(1,new_deps, ...
        matlab.codetools.requiredFilesAndProducts(p,'toponly')');
    end
    % ignore anything in matlab folder
    new_deps = new_deps(cellfun(@isempty,strfind(new_deps,matlabroot)));
    Q = setdiff(new_deps,C);
    if ~isempty(Q)
      C = cat(1,C,Q);
    end
    level= level + 1;
  end
  
  warning(w);
  
end
