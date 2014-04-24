function imwrite_gif(S,filename,varargin)
  % IMWRITE_GIF  Bootstrap `imwrite` to write an animated gif to a file given a
  % sequence of images.
  %
  % Inputs:
  %   S  h by w by c by #frames sequence of color images
  %   filename  path to .gif file
  %
  % See also: imwrite
  %

  for f = 1:size(S,4)
    [SIf,cm] = rgb2ind(S(:,:,:,f),256);
    if f==1
      imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
    else
      imwrite(SIf,cm,filename,'WriteMode','append','Delay',0);
    end
  end
end
