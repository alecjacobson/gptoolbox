function writeSVG_cubics(filename,P,C,varargin)
% writeSVG_cubics(filename,P,C,'ParameterName',ParameterValue, ...)
%
% Inputs:
%   filename  path to .svg file
%   P  #P by 2 list of control point locations
%   C  #C by 4 list of indices into P of cubics
%   Optional:
%     'MergePaths' followed by whether to merge curves into paths, default is
%       one <path> per row in C {false}
%     'Units' followed by units string {'mm'}
%   

  function append_path(fh,d_str)
      fprintf(fh,'<path style="fill:none;stroke:#000000;stroke-width:%f;stroke-miterlimit:10;" \nd="%s"></path>\n',w,d_str);
  end

  function y = range(x)
    y = max(x)-min(x);
  end

  merge_paths = false;
  units = 'mm';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MergePaths','Units'}, ...
    {'merge_paths','units'});
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

  fh = fopen(filename,'w');
  fprintf(fh,['<svg version="1.1" id="Layer_1" ' ...
    ' xmlns="http://www.w3.org/2000/svg" ' ...
    ' xmlns:xlink="http://www.w3.org/1999/xlink" ' ...
    ' width="%d%s" height="%d%s" ' ...
    ' viewBox="%d %d %d %d" ' ...
    ' xml:space="preserve">\n'], range(P(:,1)),units,range(P(:,2)),units,min(P),range(P));
  w = min(0.001*normrow(max(P)-min(P)),1);
  if merge_paths
    E = [C(:,1) C(:,end)];
    [K,A] = manifold_patches(E);
    % O(#components)
    for k = 1:max(K)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % factor out this O(#comps #curves) by pre-sorting etc.
      keep = K==k;
      Ek = E(keep,:);
      Ck = C(keep,:);
      % find path along segements and reorient if necessary
      [I,J,F] = edges_to_path(Ek);
      Ck = Ck(J,:);
      % Reorient
      Ck(F==2,:) = fliplr(Ck(F==2,:));
      assert(all(Ck(1:end-1,4) == Ck(2:end,1)))
      is_loop = I(1) == I(end);
      d_str = [sprintf('M%f,%f',P(Ck(1,1),:)')  ...
               sprintf('C%f,%f,%f,%f,%f,%f',P(Ck(1:end,2:4)',:)')];
      append_path(fh,d_str);
    end

  else
    for c = 1:size(C,1)
      d_str = sprintf('M%f,%fC%f,%f,%f,%f,%f,%f',P(C(c,:),:)');
      append_path(fh,d_str);
    end
  end
  fprintf(fh,'</svg>\n');
  fclose(fh);
end
