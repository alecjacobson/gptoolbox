function varargout = plot_cubic(C,pe,p)
  % PLOT_CUBIC  Plot a cubic Bezier curve
  %
  % [pe,p] = plot_cubic(C)
  %
  % Inputs:
  %   C  4 by dim list of control points
  % Outputs:
  %   pe  plot handles for UI
  %   p  plot handle for curve
  %
  % Example:
  %   clf;hold on;arrayfun(@(c) set(plot_cubic(P(C(c,:),:)),'Color','b'),1:size(C,1));hold off;
  t = linspace(0,1)';
  P = cubic_eval(C,t);
  if ~exist('p','var') || isempty(p)
    p = plt(P,'-k','LineWidth',2);
  else
    set(p,'XData',P(:,1),'YData',P(:,2));
  end
  ish = ishold;
  hold on;
  aiblue = hex2dec(['4D';'80';'FF'])'/255;
  if ~exist('pe','var') || isempty(pe)
    pe = plot_edges(C,[1 2;3 4],'-o','Color',aiblue,'LineWidth',1);
  else
    set(pe(1),'XData',C(1:2,1),'YData',C(1:2,2));
    set(pe(2),'XData',C(3:4,1),'YData',C(3:4,2));
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
