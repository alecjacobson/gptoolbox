demo_id = 'shape-aware'
type = 'rigid';

switch demo_id
case 'shape-aware'

  % input mesh source: *.obj, *.off
  mesh_source = 'woody.obj';
  [V,F] = load_mesh(mesh_source);
  
  % display mesh
  tsurf(F,V)
  % User clicks many times on mesh at locations of control points
  [Cx,Cy] = getpts;
  % store control points in single #P by 2 list of points
  C = [Cx,Cy];

  V = V(:,1:2);
  W = dijkstra_shepard(V,F,C);
  mlsd = MLSD2DpointsPrecompute(C',V',type,2,W');
  deform(V,F,C,'DeformMethod','mls',mlsd,'ColorScheme','MLS');
otherwise
  [im,map,alpha] = imread('~/lbs/models/woody/woody.png');
  im = ...
    im2double(im).*repmat(im2double(alpha),[1 1 size(im,3)]) + ...
    1.*repmat(1-im2double(alpha),[1 1 size(im,3)]);
  
  h = size(im,1);
  w = size(im,2);
  [X,Y] = meshgrid(1:w,1:h);
  Z = zeros(h,w);
  % display image
  handle = imshow(im);
  view(2);
  axis equal
  fprintf( ...
    ['\nCLICK on mesh at each location where you would like to add a ' ...
    'control point.\n' ...
    'Press ENTER when finished.\n\n']);
  % User clicks many times on mesh at locations of control points
  [Cx,Cy] = getpts;
  % store control points in single #P by 2 list of points
  C = [Cx,Cy];
  
  step = 10;
  [X,Y] = meshgrid(1:step:size(im,2),1:step:size(im,1));
  V = [X(:), Y(:)];
  type = 'rigid';
  
  mlsd = MLSD2DpointsPrecompute(C',V',type,2);
  deform([],[],C,'DeformMethod','mls',mlsd,'Image',im,X,Y,'ColorScheme','MLS');
end

