function [V,F] = bwstencil(im,varargin)
  % BWSTENCIL Create a 3D-printable stencil from an input binary image.
  %
  % [V,F] = bwstencil(im)
  % [V,F] = bwstencil(im,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   im  h by w by 1 binary image
  %   Optional inputs (see bwmesh.m)
  % Outputs:
  %   V  #V by 3 list of vertex positions with V(:,3) = 0 or 1 for top or
  %   bottom layer
  % 
  % Example:
  %   % Read in image, convert to binary, shrink
  %   im = imresize(im2bw(imread('hans-hass.jpg')),0.25);
  %   % Pad image to create boundary frame
  %   im = padarray(im,ceil(max(size(im))*0.125*[1;1]));
  %   % Generate stencil of thickness = 1px
  %   [V,F] = bwstencil(im,'Tol',1.5,'SmoothingIters',50);
  %   % Change thickness to 10% size of image
  %   V(:,3) = V(:,3) * 0.1*max(size(im));
  %   % render result:
  %   subplot(2,1,1);
  %   imshow(im);
  %   subplot(2,1,2);
  %   tsurf(F,V,'EdgeColor','none','FaceColor',[0.3 0.4 1.0],'FaceAlpha',0.9);
  %   l = light('Style','infinite','Position',[0.3 0.2 1]);
  %   axis equal;
  %
  % See also: bwmesh
  % 
  if ~islogical(im)
    warning('Converting non-binary input to binary');
    im = im2bw(im);
  end
  imh = imfill(im,'holes');
  if any(im(:) ~= imh(:))
    warning('input contained islands, removing...');
    im = imh;
  end
  im = imresize(im,2,'nearest');
  [V,F] = bwmesh(~im,varargin{:});
  [V,F] = extrude(V,F);
  [V,F] = remesh_planar_patches(V,F);
end
