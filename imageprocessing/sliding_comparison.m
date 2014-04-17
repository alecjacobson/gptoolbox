function [ph,ih] = sliding_comparison(A,B)
  % SLIDING_COMPARISON Create a figure which allows the user to mouse over and
  % move a sliding splitter to compare two images interactively.
  %
  % Inputs:
  %   A  h by w by c left image (should be double)
  %   B  h by w by c right image (should be double)
  % Outputs:
  %   ph  handle to plot of separator line
  %   ih  handle to image line
  %
  function C = split(A,B,x)
    x = max(min(x,size(A,2)),1);
    C = [A(:,1:round(x)-1,:) B(:,round(x):end,:)];
  end

  function onmove(src,ev)
    % get current mouse position, and remember old one
    pos=get(gca,'currentpoint');
    x = pos(1,1,1);
    if ishandle(ih) && ishandle(ph)
      set(ph,'XData',repmat(x,1,2));
      set(ih,'CData',split(A,B,x));
    else
      % Clean up
      set(gcf,'windowbuttonmotionfcn',old_onmove);
    end
    drawnow
  end

  w = size(A,2);
  h = size(A,1);
  % `imshow` will clamp but `set(...,'CData',...)` will not
  clamp = @(X) max(min(X,1),0);
  A = clamp(A);
  B = clamp(B);

  x = w/2;
  ih = imshow(split(A,B,x));
  hold on;
    ph = plot([x;x],[0;h],'LineWidth',2);
  hold off;

  old_onmove = get(gcf,'windowbuttonmotionfcn');
  set(gcf,'windowbuttonmotionfcn',@onmove);
end
