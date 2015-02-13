function f = figgif(filename)
  frame = getframe(gcf);
  [SIf,cm] = rgb2ind(frame.cdata,256);
  if ~exist(filename,'file')
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
  end
end
