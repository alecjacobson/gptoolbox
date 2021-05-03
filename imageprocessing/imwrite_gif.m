function imwrite_gif(S,filename,varargin)
  % IMWRITE_GIF  Bootstrap `imwrite` to write an animated gif to a file given a
  % sequence of images.
  %
  % Inputs:
  %   S  h by w by c by #frames sequence of color images
  %   filename  path to .gif file
  %
  % Example:
  %   % Read an animation
  %   [X,M] = imread('peaks.gif');
  %   % Convert to raw color image sequence
  %   Y = permute(reshape(M(X(:)+1,:),[size(X,1),size(X,2),size(X,4),3]),[1 2 4 3]);
  %   % Trim animation to fit
  %   C = imtrim(Y);
  %   % Write back to animated .gif
  %   imwrite_gif(C(:,:,:,[1:end end-1:-1:1]),'peaks2.gif');
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
