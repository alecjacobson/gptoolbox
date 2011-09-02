function patcht(FF,VV,TF,VT,I,Options)
% This function PATCHT, will show a triangulated mesh like Matlab function
% Patch but then with a texture.
%
% Alec: Essentially this calls "surface" for each triangle of the triangle mesh
% with the appropriate texture coordinates
%
% patcht(FF,VV,TF,VT,I,Options);
%
% inputs,
%   FF : Face list #F x 3 with vertex indices
%   VV : Vertices #V x 3
%   TF : Texture list #F x 3 with texture vertex indices
%   VT : Texture Coordinates s 2 x #VT, range must be [0..1] or real pixel postions
%   I : The texture-image RGB [O x P x 3] or Grayscale [O x P] 
%   Options : Structure with options for the textured patch such as
%           EdgeColor, EdgeAlpha see help "Surface Properties :: Functions"
%
%   Options.PSize : Special option, defines the image texturesize for each 
%           individual  polygon, a low number gives a more block 
%           like texture, defaults to 64;
%
% note: 
%   On a normal PC displaying 10,000 faces will take about 6 sec.
%
% Example,
%
%  % Load Data;
%   load testdata;
%  % Show the textured patch
%   figure, patcht(FF,VV,TF,VT,I);
%  % Allow Camera Control (with left, right and center mouse button)
%   mouse3d
%
%  % file containing image
%  image_source = 'woody.png';
%  % load image with alpha mask
%  [im,map,alpha] = imread(image_source);
%  % mesh based on alpha values with approximately 100 boundary vertices
%  [V,F] = png2mesh(image_source,0,100);
%  % treat alpha as white in display
%  im = over(im,alpha,ones(size(im)),1);
%  % dummy display white image the same size as input image
%  imshow(ones([max(V(:,2))-min(V(:,2)) max(V(:,1))-min(V(:,1))]));
%  % display original mesh
%  % reverse y-direction so mesh display is upright (we're displaying the mesh
%  % ontop of the imshow() display to get pixel perfect images)
%  patcht(F,[V(:,1) size(im,1)-V(:,2)],F,V,permute(flipdim(im,1),[2 1 3]));
%  % create some trivial deformation
%  new_V = [V(:,1) 2*V(:,2)];
%  figure;
%  % dummy display white image the same size as input image
%  imshow(ones([max(new_V(:,2))-min(new_V(:,2)) max(new_V(:,1))-min(new_V(:,1))]));
%  % display deformed mesh
%  patcht(F,[new_V(:,1) max(new_V(:,2)+0.5)-new_V(:,2)],F,V,permute(flipdim(im,1),[2 1 3]));
%
% Function is written by D.Kroon University of Twente (July 2010)

% FaceColor is a texture
Options.FaceColor='texturemap';
% noe edges
Options.EdgeColor='none';

% Size of texture image used for every triangle
if(isfield(Options,'PSize'))
    sizep=round(Options.PSize(1));
    Options=rmfield(Options,'PSize');
else
    sizep=64;
end

% Check input sizes
if(size(FF,2)~=size(TF,2))
    error('patcht:inputs','Face list must be equal in size to texture-index list');
end

if((ndims(I)~=2)&&(ndims(I)~=3))
    error('patcht:inputs','No valid Input texture image');
end

% Detect if grayscale or color image
switch(size(I,3))
    case 1
        iscolor=false;
    case 3
        iscolor=true;
    otherwise
        error('patcht:inputs','No valid Input texture image');
end

   
if(max(VT(:))<2)
    % Remap texture coordinates to image coordinates
    VT2(:,1)=(size(I,1)-1)*(VT(:,1))+1;
    VT2(:,2)=(size(I,2)-1)*(VT(:,2))+1;
else
    VT2=VT;
end

% Calculate the texture interpolation values
[lambda1 lambda2 lambda3 jind]=calculateBarycentricInterpolationValues(sizep);
 
% Split texture-image in r,g,b to allow fast 1D index 
Ir=I(:,:,1); if(iscolor), Ig=I(:,:,2); Ib=I(:,:,3); end

% The Patch used for every triangle (rgb)
Jr=zeros([(sizep+1) (sizep+1) 1],class(I));
if(iscolor)
    Jg=zeros([(sizep+1) (sizep+1) 1],class(I));
    Jb=zeros([(sizep+1) (sizep+1) 1],class(I));
end

% if no z coordinates are given then use 0s
if size(VV,2) < 3
  VV = [VV zeros(size(VV,1),1)];
end

hold on;
% Loop through all triangles of the mesh
for i=1:size(FF,1)
    % Get current triangle vertices and current texture-vertices
    V=VV(FF(i,:),:); 
    Vt=VT2(TF(i,:),:); 
    
    % Define the triangle as a surface
    x=[V(1,1) V(2,1); V(3,1) V(3,1)];
    y=[V(1,2) V(2,2); V(3,2) V(3,2)];
    z=[V(1,3) V(2,3); V(3,3) V(3,3)];
    
    % Define the texture coordinates of the surface
    tx=[Vt(1,1) Vt(2,1) Vt(3,1) Vt(3,1)];
    ty=[Vt(1,2) Vt(2,2) Vt(3,2) Vt(3,2)] ;
    xy=[tx(1) ty(1); tx(2) ty(2); tx(3) ty(3); tx(3) ty(3)];

    % Calculate texture interpolation coordinates
    pos(:,1)=xy(1,1)*lambda1+xy(2,1)*lambda2+xy(3,1)*lambda3;
    pos(:,2)=xy(1,2)*lambda1+xy(2,2)*lambda2+xy(3,2)*lambda3;
    pos=round(pos); pos=max(pos,1); pos(:,1)=min(pos(:,1),size(I,1)); pos(:,2)=min(pos(:,2),size(I,2));
    posind=(pos(:,1)-1)+(pos(:,2)-1)*size(I,1)+1;
    
    % Map texture to surface image
    Jr(jind)=Ir(posind);
    J(:,:,1)=Jr; 
    if(iscolor)
        Jg(jind)=Ig(posind); 
        Jb(jind)=Ib(posind);
        J(:,:,2)=Jg; 
        J(:,:,3)=Jb;
    end
    
    % Show the surface
    surface(x,y,z,J,Options)
end
hold off; 

function [lambda1 lambda2 lambda3 jind]=calculateBarycentricInterpolationValues(sizep)
% Define a triangle in the upperpart of an square, because only that
% part is used by the surface function
x1=sizep; y1=sizep; x2=sizep; y2=0; x3=0 ;y3=0;
% Calculate the bary centric coordinates (instead of creating a 2D image
% with the interpolation values, we map them directly to an 1D vector)
detT = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
[x,y]=ndgrid(0:sizep,0:sizep); x=x(:); y=y(:);
lambda1=((y2-y3).*(x-x3)+(x3-x2).*(y-y3))/detT;
lambda2=((y3-y1).*(x-x3)+(x1-x3).*(y-y3))/detT;
lambda3=1-lambda1-lambda2;
% Make from 2D (surface)image indices 1D image indices
[jx jy]=ndgrid(sizep-(0:sizep)+1,sizep-(0:sizep)+1);
jind=(jx(:)-1)+(jy(:)-1)*(sizep+1)+1;



    
 
