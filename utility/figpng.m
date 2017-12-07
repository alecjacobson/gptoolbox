function f = figpng(filename,varargin)
  % FIGPNG Write the current figure to filename. If filename exists, replace it
  %
  % Inputs:
  %   filename  path to png file
  % Outputs:
  %   f  flag whether file already existed
  % 
  frame = getframe(gcf);
  f = exist(filename,'file');
  imwrite(frame.cdata,filename);
end

