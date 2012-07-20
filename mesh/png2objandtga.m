function [V,F,img] = png2objandtga( ...
  filename, ...
  laplacian_smoothness_iterations, ...
  max_points_on_boundary)
% PNG2OBJANDTGA triangulate a png file based on its alpha mask and spit out
% corresponding tga (texture) and obj (mesh) files
%
% png2objandtga(filename)
% png2objandtga(filename,laplacian_smoothness_iterations,max_points_on_boundary)
%
% Inputs:
%   filename  path to .png file
%   Optional:
%     laplacian_smoothness_iterations  number of laplacian smoothness
%       iterations to perform on boundary before meshing {0}
%     max_points_on_boundary {0.5(#width+#height)}
% Outputs:
%   V  #V by 2 list of mesh vertex positions
%   F  #F by 3 list of triangle indices
%   img  image
%


  basename = regexprep(filename,'\.png$','');
  img = imread(filename);
  if(0 ~= mod(size(img,1),4) || 0 ~= mod(size(img,2),4))
    warning(['Image dimensions are not a multiple of 4...' ...
      ' tga may not load in openGL properly']);
  end
  
  
  
if ~exist('laplacian_smoothness_iterations','var')
  laplacian_smoothness_iterations = 0;
end

if ~exist('max_points_on_boundary','var')
  max_points_on_boundary = 0.5*sum(size(img));
end

  [V,F] = png2mesh( ...
    filename, ...
    laplacian_smoothness_iterations, ...
    max_points_on_boundary);
  tsurf(F,V);
  V = [V(:,1) V(:,2) 0*V(:,1)];
  writeOBJ([basename '.obj'],V,F);
  % convert png to tga
  path_to_convert = '/opt/local/bin/convert';
  [status, result] = ...
    system(['export DYLD_LIBRARY_PATH="";' path_to_convert ' -flatten ' basename '.png' ' ' basename '.tga']);
  [status, result] = ...
    system(['export DYLD_LIBRARY_PATH="";' path_to_convert ' +matte ' basename '.tga' ' ' basename '.tga']);
end
