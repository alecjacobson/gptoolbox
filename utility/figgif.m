function f = figgif(filename,varargin)
  % FIGGIF Write the current figure to filename. If filename exists, append it
  % as an animation frame with zero delay
  %
  % Inputs:
  %   filename  path to gif file
  % Outputs:
  %   f  flag whether file already existed
  % 
  frame = getframe(gcf);
  [SIf,cm] = rgb2ind(frame.cdata,256,varargin{:});
  f = exist(filename,'file');
  if ~f
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
  end
end
