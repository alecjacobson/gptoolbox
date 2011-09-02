function interactive_fit_plot(x,y,x_min,x_max,y_min,y_max,x_label,y_label,log_base,m,b)
  hold off;
  % scatter plot
  plot(x,y,'o','MarkerEdgeColor',[0.2 0.5 0.2],'MarkerSize',10,'LineWidth',1.5);%'MarkerFaceColor',[0.2 0.5 0.2 ]);
  axis([x_min,x_max,y_min,y_max]);
  set(gca,'FontSize',30);
  % hack to get yticks off of xticks
  yticks = y_min:y_max;
  yticks = yticks+(yticks(2)-yticks(1))/16;
  % further hack to move rightmost xtick to left of right outline
  xticks = x_min:x_max;
  xticks(2) = xticks(2)-(xticks(2)-xticks(1))/8;
  set(gca,'YTick',yticks);
  set(gca,'XTick',xticks);

  % line
  if(exist('m') && exist('b'))
    fit = m*x+b;
    hold on;
    plot(x,fit,'LineWidth',2,'Color','k');
    title(['Slope: ' num2str(m) ]);
  end

  if(exist('x_label'))
    %xlabel(x_label);
  end

  if(exist('y_label'))
    %ylabel(y_label);
  end


  if(exist('log_base') & log_base)
    set(gca,'XTickLabel', ...
      [repmat('    1.E',size((x_min:x_max),2),1) num2str((x_min:x_max)') ...
      ]);
    set(gca,'YTickLabel', ...
      [repmat('1.E',size((y_min:y_max),2),1) num2str((y_min:y_max)') ...
      ]);
  end
end
