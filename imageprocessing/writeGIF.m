function writeGIF(filename,X,M,D)
  % WRITEGIF imwrite is bad at writing animated figs. This allows per frame
  % delays
  % 
  % writeGIF(filename,X,M,D)
  %
  % Inputs:
  %   filename  path to .gif file
  %   X  h by w by 3 by #frames list of indexed images
  %   M  h by w by (1|#frames) list of maps
  %   D  (1|#frames)  list of delays
  % 
  % See also, hack to get imread to read animated gifs correctly:
  % http://www.alecjacobson.com/weblog/?p=4726
  %
  % Example:
  %  % Load in a one-way animation
  %  [X,M] = imread('bunny.gif');
  %  % Trim to animations extent
  %  Y = imtrim(X,'Map',M,'Threshold',0);
  %  % Write a "boomerang": looping back-and-forth image, pausing at endpoints
  %  % for 0.5 seconds
  %  writeGIF( ...
  %    'output.gif', ...
  %    Y(:,:,:,[1:end end end:-1:1]), ...
  %    M(:,:,[1:end end end:-1:1]), ...
  %    sparse([size(Y,4)+1 2*size(Y,4)+1],1,0.5,2*size(Y,4)+1,1)); 
  % 
  imwrite(X(:,:,:,1),M(:,:,1),filename,'DelayTime',D(1));
  arrayfun(@(f) imwrite(X(:,:,:,f),M(:,:,min(f,end)),filename,'Delay',D(min(f,end)),'WriteMode','append'),2:size(X,4));
end
