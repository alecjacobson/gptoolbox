function p = add_isolines(T,varargin);
  % ADD_ISOLINES Add isolines for plotted mesh (t) with scalar per-vertex
  % pseudocolor data
  %
  % p = add_isolines()  % Apply to all 
  % h = add_shadow(T)
  % h = add_shadow(T,ParametersPassedToPlot, ...)
  %
  % InpuT:
  %   T  #T list of trisurf handles {[] --> find all trisurf in `gca`}
  % OutpuT:
  %   p  #T list of isolines plot handles
  %
  if ~exist('T','var') || isempty(T)
    c = get(gca,'Children');
    T = c(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Patch'),c));
  end
  if iscell(T)
    T = [T{:}];
  end
  p = {};
  for ti = 1:numel(T)
    t = T(ti);
    n = size(colormap,1)+1;
    a = caxis;
    iso = linspace(a(1),a(2),n);
    V = t.Vertices;
    F = t.Faces;
    S = t.FaceVertexCData;
    assert(size(S,2) == 1);
    [LS,LD,I] = isolines(V,F,S,iso);
    hold on;
    p{end+1} = plot3([LS(:,1) LD(:,1)]',[LS(:,2) LD(:,2)]',[LS(:,3) LD(:,3)]','Color','k',varargin{:});
    hold off;
  end
end
