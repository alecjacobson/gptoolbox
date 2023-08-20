function [W,D,E] = cubic_winding_number(C,V,ispoly)
  % CUBIC_WINDING_NUMBER Cubic the winding number of points in V with respect to
  % a cubic Bezier curve in C
  %
  % [W,D,E] = cubic_winding_number(C,V)
  %
  % Inputs:
  %   C  4 by 2 list of control points
  %   V  #V by 2 list of query points
  % Outputs:
  %   W  #V list of winding number values
  %   D  #V list of max-depths of recursive algorithm 
  %   E  #V list of total evaluations of recursive algorithm
  % 

  function I = inpolygon_convex(V,C);
    %I = inpolygon(V(:,1),V(:,2),C(:,1),C(:,2));
    %I = abs(winding_number(C,[1:size(C,1);2:size(C,1) 1]',V))>0.5;
    %% Geez, this is barely faster than winding_number
    %Nx = C([2:end 1],2)-C(:,2);
    %Ny = C(:,1)-C([2:end 1],1);
    %S = (V(:,1)-C(:,1)').*Nx' + (V(:,2)-C(:,2)').*Ny';
    %I = all(S<0,2);
    % Oh, matlab, why:
    I = all(((V(:,1)-C(:,1)').*(C([2:end 1],2)-C(:,2))' + (V(:,2)-C(:,2)').*(C(:,1)-C([2:end 1],1))')<0,2);
  end

  if nargin<3 || isempty(ispoly)
    ispoly = false;
  end

  max_depth = 40;
  W = zeros(size(V,1),1);
  D = zeros(size(V,1),1);
  E = zeros(size(V,1),1);
  k = 0;
  % "queue"
  Q = {{C,(1:size(V,1))',1}};
  while ~isempty(Q)
    % Pop off back: faster but uglier traversal for animation
    Qi = Q{end};
    Ci = Qi{1};
    Ji = Qi{2};
    depth = Qi{3};
    Q = Q(1:end-1);
    %% Pop off front: better looking traversal for animation, but slower
    %Qi = Q{1};
    %Ci = Qi{1};
    %Ji = Qi{2};
    %depth = Qi{3};
    %Q = Q(2:end);
    D(Ji) = max(D(Ji),depth);
    E(Ji) = E(Ji)+1;
    if depth >= max_depth
      %warning('spline_winding_number exceeded max depth (%d)\n',max_depth);
      continue;
    end
    if size(Ci,1)<3
      I = false(numel(Ji),1);
    else
      H = convhull(Ci);
      I = inpolygon_convex(V(Ji,:),Ci(H(1:end-1),:));
    end
    Wi = winding_number(Ci,[size(Ci,1) 1],V(Ji(~I),:));
    W(Ji(~I)) = W(Ji(~I)) - Wi;
    %WW = zeros(size(W));
    %WW(Ji(~I)) = -Wi;
    %ss = [915 938];
    %surf( ...
    %  reshape(V(:,1),ss), ...
    %  reshape(V(:,2),ss), ...
    %  reshape(0*V(:,2),ss), ...
    %  'CData', ...
    %  reshape(WW,ss), fphong, 'EdgeColor', 'none');
    %hold on
    %[pe,p] = plot_cubic(C);
    %set(p,'Color',0.75*[1 1 1]);
    %arrayfun(@(p) set(p,'Color',1+0.5*(get(p,'Color')-1)),pe);
    %plot_cubic(Ci);
    %plot_edges(Ci,[1 4],':','LineWidth',2,'Color',0.4*[1 1 1]);
    %hold off;
    %view(2);
    %axis equal;
    %axis tight;
    %caxis([-1 1])
    %colormap((flipud(cbrewer('RdBu',45))))
    %set(gca,'YDir','reverse');
    %set(gcf,'color','w');
    %set(gca,'Visible','off','Position',[0 0 1 1]);
    %figpng(sprintf('cubic-winding_number-%02d.png',k));
    %k=k+1;

    Ki = Ji(I);
    if ~isempty(Ki)
      if ispoly
        s = ceil(size(Ci,1)/2);
        C1 = Ci(1:s,:);
        C2 = Ci(s:end,:);
        Q{end+1} = {C1,Ki,Qi{3}+1};
        Q{end+1} = {C2,Ki,Qi{3}+1};
      else
        [C1,C2] = cubic_split(Ci,0.5);
        Q{end+1} = {C1,Ki,Qi{3}+1};
        Q{end+1} = {C2,Ki,Qi{3}+1};
      end
    end
  end
end
