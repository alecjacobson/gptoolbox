function [PP,CC,I] = spline_upsample(P,C,varargin)
  % SPLINE_UPSAMPLE Upsample a spline 
  %
  % [PP,CC,I] = spline_upsample(P,C);
  %
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  %   Optional:
  %     'OnlySelected' see `upsample.m`
  %     'Iterations' followed by number of iterations to conduct {1}
  % Outputs:
  %   PP  #PP by dim list of control point locations
  %   CC  #CC by 4 list of indices into PP of cubic Bezier curves
  %   I  #CC list of indices into (1:size(C,1)) revealing birth segments
  %
  % Example:
  %   [uP,uC] = spline_upsample(dP,dC,'Iterations',inf,'OnlySelected',@(P,C) edge_lengths(P,C(:,1:2))+edge_lengths(P,C(:,2:3))+edge_lengths(P,C(:,3:4)) > 50);
  iters = 1;
  sel = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'OnlySelected','Iterations'}, ...
    {'sel','iters'});
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

  if isempty(sel)
    sel = (1:size(C,1))';
  end

  if iters < 1
    PP = P;
    CC = C;
    I = (1:size(C,1))';
    return;
  end

  sel_fun = [];
  if isa(sel,'function_handle')
    sel_fun = sel;
    sel = sel_fun(P,C);
  end
  if ~islogical(sel)
    sel = accumarray(sel,true,[size(C,1),1]);
  end
  sel = sel(:);
  if ~any(sel)
    PP = P;
    CC = C;
    I = (1:size(C,1))';
    return;
  end

  PP = [];
  CC = [];
  I = [];
  for c = 1:size(C,1)
    if sel(c)
      [Pc,Cc] = cubic_subdivide(P(C(c,:),:),0.5);
    else
      % This is a bit silly with a lot of copying. And makes this not output
      % sensitive.
      Pc = P(C(c,:),:);
      Cc = [1 2 3 4];
    end
    I = [I;repmat(c,size(Cc,1),1)];
    CC = [CC size(PP,2)+Cc'];
    PP = [PP Pc'];
  end
  PP = PP';
  CC = CC';
  [PP,~,~,CC] = remove_duplicate_vertices(PP,0,'F',CC);

  if iters == 1
    return;
  end

  [PP,CC,Ir] = spline_upsample(PP,CC,varargin{:},'Iterations',iters-1);
  I = I(Ir);


end
