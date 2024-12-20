function varargout = plot_spline(P,C,varargin)
  % PLOT_SPLINE Plot a cubic Bezier spline.
  %
  % [pe,p] = plot_spline(P,C)
  %
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  % Outputs:
  %   pe  plot handles for UI
  %   p  plot handle for curves
  % Example:
  %   cellfun(@(pe) arrayfun(@(p) set(p,'Color','r'),pe),plot_spline(P,C))

  if numel(varargin) == 1 && iscell(varargin{1}) && numel(varargin{1}) == 2
    p_varargin = varargin{1}{2};
    pe_varargin = varargin{1}{1};
  else
    p_varargin = {};
    pe_varargin = varargin;
  end
  
  assert(max(C(:))<=size(P,1));
  assert(min(C(:))>=1);
  p = {};
  pe = {};
  ish = ishold;
  for c = 1:size(C,1)
    [pe{c},p{c}] = plot_cubic(P(C(c,:),:),[],[],{pe_varargin,p_varargin});
    hold on;
  end
  hold off;
  if ish
    hold on
  end
  if nargout>=1
    varargout{1} = pe;
    if nargout >= 2
      varargout{2} = p;
    end
  end
end

