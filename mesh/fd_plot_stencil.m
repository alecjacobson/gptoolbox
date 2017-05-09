function [s,t] = fd_plot_stencil(side,V)
  % Visualization of the stencil (non-zeros in a row V of a matrix).
  %
  % [s,t] = fd_plot_stencil(side,V)
  %
  % Inputs:
  %   side  [h w]
  %   V  h*w by 1 row of a sparse matrix
  % Outputs:
  %   s  handle to surf plot
  %   t  handle to text plot
  %
  % Example:
  %   Q = fd_laplacian([5 5]);
  %   fd_plot_stencil([5 5],Q(ceil(end/2),:));
  %
  h = side(1);
  w = side(2);
  [X,Y] = meshgrid((1:w+1)-0.5,(1:h+1)-0.5);
  s = surf(X,Y,0*X,reshape(V,side),'EdgeColor','k');
  [X,Y] = meshgrid((1:w),(1:h));
  t = text(X(:),Y(:),num2str(V(:),'%0.1f'), ...
      'HorizontalAlignment','center', ...
      'FontSize',20, ...
      'BackgroundColor',[0.7 0.7 0.7]);
  view(2);
  axis equal;
end

