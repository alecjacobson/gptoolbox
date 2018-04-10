function [h,L,M,ground] = add_shadow(T,L,varargin)
  % ADD_SHADOW  Add a shadow for plotted mesh (t) according to light (l)
  %
  % h = add_shadow()  % Apply to all 
  % h = add_shadow(T,L)
  % h = add_shadow(T,L,'ParameterName',ParameterValue,...)
  % h = add_shadow([],[],'ParameterName',ParameterValue,...) % apply to all
  %
  % Inputs:
  %   T  #T list of trisurf handles {[] --> find all trisurf in `gca`}
  %   L  #L list of lights {[] --> find all light in `gca`}
  %   Optional:
  %     'Ground'  ground plane equation {[0 0 -1 min(Z)]}
  %     'Nudge'  nudge the ground plane down a bit
  %     'Color' followed by 3-vector color {get(gcf,'Color')*0.9}
  %     'BackgroundColor' followed by 3-vector color {get(gcf,'Color')*0.9}
  %     'Fade'  followed by:
  %        'none' constant shadow color
  %        'local' fade darker away from contact with ground (ape a spotlight)
  %        {'infinite'} fade lighter away from contact ground (ape infinite
  %          light)
  % Outputs:
  %   h  #T*#L list of output shadow trisurf handles
  %   L  #L list of lights
  %   M  4 by 4 by #T*#L shadow projection matrices
  %
  % Example:
  %   t = tsurf(F,V,'EdgeColor','none',fphong,'SpecularStrength',0.1);
  %   l = [ 
  %      light('Position',[-10 -10 13],'Style','local');
  %      light('Position',[10 -10  13],'Style','local')];
  %   camproj('persp');
  %   axis equal;
  %   h = add_shadow(t,l);
  %   apply_ambient_occlusion(t);
  %

  % default values
  ground = [];
  nudge = 0;
  color = get(gcf,'Color')*0.9;
  background_color = get(gcf,'Color');
  fade = 'infinite';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Ground','Nudge','BackgroundColor','Color','Fade'}, ...
    {'ground','nudge','background_color','color','fade'});
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

  if ~exist('T','var') || isempty(T)
    c = get(gca,'Children');
    T = c(arrayfun(@(x) ...
      isa(x,'matlab.graphics.primitive.Patch') || ...
      isa(x,'matlab.graphics.chart.primitive.Scatter') ...
      ,c));
  end
  if iscell(T)
    T = [T{:}]';
  end
  T = T(:);

  if ~exist('L','var') || isempty(L)
    c = get(gca,'Children');
    L = c(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Light'),c));
    if isempty(L)
      L = camlight;
    end
  end

  if isempty(ground)
    minZ = inf;
    for t = T'
      switch class(t)
      case 'matlab.graphics.primitive.Patch'
        V = t.Vertices;
      case 'matlab.graphics.chart.primitive.Scatter'
        V = [t.XData;t.YData;t.ZData]';
      otherwise
        class(t)
      end
      minZ = min([V(:,3);minZ]);
    end
    ground = [0 0 -1 minZ];
  end
  ground(4) = ground(4)-nudge;

  h = {};
  % need to specify that there are 0 "tubes"
  M = zeros([0,0,0]);
  for t = T'
    switch class(t)
    case 'matlab.graphics.primitive.Patch'
      V = t.Vertices;
    case 'matlab.graphics.chart.primitive.Scatter'
      V = [t.XData;t.YData;t.ZData]';
    otherwise
      class(t)
    end
    for l = L'
      % plane equation
      % 0 = ax + by + cz + d
      light_pos = [l.Position strcmp(l.Style,'local')];
      d = ground * light_pos';
      shadow_mat = d*eye(4) - light_pos'*ground;
      U = [V ones(size(V,1),1)]*shadow_mat';
      U = bsxfun(@rdivide,U(:,1:3),U(:,4));
    
      hold on;
      switch class(t)
      case 'matlab.graphics.primitive.Patch'
        tsh = trisurf(t.Faces,U(:,1),U(:,2),U(:,3), ...
          'FaceColor',color, ...
          'DiffuseStrength',0,'SpecularStrength',0, ...
          'AmbientStrength',1, ...
          'EdgeColor','none');
      case 'matlab.graphics.chart.primitive.Scatter'
        tsh = copyobj(t,t.Parent);
        tsh.XData = U(:,1);
        tsh.YData = U(:,2);
        tsh.ZData = U(:,3);
        tsh.MarkerFaceColor = color;
        tsh.MarkerEdgeColor = color;
      end
      hold off;
      switch fade
      case {'local','infinite'}
        D = matrixnormalize( ...
          sum(bsxfun(@times,U(:,1:2),l.Position(1:2)),2));
        switch fade
        case 'infinite'
          D = 1.0-D;
        end
        ca = caxis;
        C = bsxfun(@plus,color,bsxfun(@times,D,background_color-color));
        switch class(t)
        case 'matlab.graphics.primitive.Patch'
          tsh.FaceVertexCData = C;
          tsh.FaceColor = 'interp';
        case 'matlab.graphics.chart.primitive.Scatter'
          tsh.MarkerEdgeColor = 'flat';
          tsh.MarkerFaceColor = 'flat';
          tsh.CData = C;
        end
        caxis(ca);
      end
      h = {h{:} tsh};
      M(:,:,end+1) = shadow_mat;
    end
  end
end
