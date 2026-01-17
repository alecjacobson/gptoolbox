function [Lines,Vertices,Objects]=isocontour(I,isovalue)
% This function ISOCONTOUR computes the isocontour geometry for 
% a certain 2D image and isovalue. To Extract the isocontour geometry it 
% uses Marching Squares and linear interpolation. Followed by sorting the
% contour geometry into separate sorted contours.
%
% This function is the 2D equivalent of Isosurface extraction 
% using Marching Cubes in 3D.
%
%
%   [Lines,Vertices,Objects]=isocontour(I,isovalue);
%
% inputs,
%   I : An 2D image (grey-scale)
%   isovalue : The Iso-value of the contour
%
% outputs,
%   Lines : An array describing all the Line-pieces of the isocontour
%           geomtery, with a N x 2 index list of vertices
%   Vertices : Vertices (Corners) of the lines M x 2 list of X,Y
%           coordinates
%   Objects : A 1 x K  cell array with in every cell a list of indices
%           corresponding to one connect isocontour. If the isocontour
%           is closed then the last index value is equal to first index
%           value.
%
% Note : This function displays the image with isocontours if no output
%       is defined.
%
% Example,
%    I = im2double(imread('rice.png'));
%    isocontour(I,0.5);
%
% Example,
%    I = im2double(imread('rice.png'));
%    [Lines,Vertices,Objects]=isocontour(I,0.5);
%    figure('renderer','opengl'), imshow(I), hold on;
%    for i=1:length(Objects)
%         Points=Objects{i};
%         plot(Vertices(Points,2),Vertices(Points,1),'Color',rand(3,1));
%    end
% 
% Example,
%    I = im2double(imread('rice.png'));
%    [Lines,Vertices]=isocontour(I,0.5);
%    figure, imshow(I), hold on;
%    V1=Vertices(Lines(:,1),:); V2=Vertices(Lines(:,2),:);
%    plot([V1(:,2) V2(:,2)]',[V1(:,1) V2(:,1)]','b');
%    % Calculate Contour Normals
%    N = V1-V2; L = sqrt(N(:,1).^2+N(:,2).^2)+eps; 
%    N(:,1)=N(:,1)./L; N(:,2)=-N(:,2)./L;
%    plot([V1(:,2) V1(:,2)+N(:,1)*5]',[V1(:,1) V1(:,1)+N(:,2)*5]','r');
%
%
% This function is written by D.Kroon University of Twente (March 2011)  


% Check Inputs
if(nargin==0), error('isocontour:input','no input image defined'); end
if(ndims(I)~=2), error('isocontour:input','image must be 2D'); end
if(size(I,3)>1), error('isocontour:input','image cannot be RGB'); end
if(nargin==1), isovalue=0.5; end

% Get the Line Pieces
[Vertices,Lines]=LookupDB(I,isovalue);

% Sort the Line Pieces into objects
if(nargout==3||nargout==0)
    Objects=SortLines2Objects(Lines);
end

% Show image, if no output asked
if(nargout==0)
    figure, imshow(I), hold on;
    for i=1:length(Objects)
         Points=Objects{i};
         plot(Vertices(Points,2),Vertices(Points,1),'Color',[1 0 0]);
    end
end


function Objects=SortLines2Objects(Lines)
% Object index list
Obj=zeros(100,2); Obj(1,1)=1; nObjects=1; reverse=false;
for i=1:size(Lines,1)-1
    F=Lines(i,2);
    Q=i+1;
    R=find(any(Lines(Q:end,:)==F,2),1,'first');
    if(~isempty(R))
        R=R+i;
        TF=Lines(Q,:);
        if(Lines(R,1)==F),Lines(Q,:)=Lines(R,[1 2]); else Lines(Q,:)=Lines(R,[2 1]); end
        if(R~=Q), Lines(R,:)=TF; end
    else
        F=Lines(Obj(nObjects,1),1);
        R=find(any(Lines(Q:end,:)==F,2),1,'first');
        if(~isempty(R))
            reverse=true;
            Lines(Obj(nObjects,1):i,:)=Lines(i:-1:Obj(nObjects,1),[2 1]);
            R=R+i;
            TF=Lines(Q,:);
            if(Lines(R,1)==F),Lines(Q,:)=Lines(R,[1 2]); else Lines(Q,:)=Lines(R,[2 1]); end
            if(R~=Q), Lines(R,:)=TF; end
        else
            if(reverse)
                Lines(Obj(nObjects,1):i,:)=Lines(i:-1:Obj(nObjects,1),[2 1]);
                reverse=false;
            end
            Obj(nObjects,2)=i;
            nObjects=nObjects+1;
            Obj(nObjects,1)=i+1;
        end
    end
end
Obj(nObjects,2)=i+1;

% Object index list, to real connect object lines
Objects=cell(1,nObjects);
for i=1:nObjects
    % Determine if the line is closed
    if(Lines(Obj(i,1),1)==Lines(Obj(i,2),2)) 
        Objects{i}=[Lines(Obj(i,1):Obj(i,2),1);Lines(Obj(i,1),1)];
    else
        Objects{i}=Lines(Obj(i,1):Obj(i,2),1);
    end
end

function [V,F]=LookupDB(Img,isovalue)

% Describe the base-polygons by edges
%  Edge number
%   1
%  ***
% 0* *2
%  ***
%   3
%
I=zeros(16,4);
I(1,:)= [4 4 4 4]; % [0 0;0 0]
I(2,:)= [1 0 4 4]; % [1 0;0 0] 
I(3,:)= [2 1 4 4]; % [0 1;0 0] 
I(4,:)= [2 0 4 4]; % [1 1;0 0] 
I(5,:)= [0 3 4 4]; % [0 0;1 0] 
I(6,:)= [1 3 4 4]; % [1 0;1 0] 
I(7,:)= [2 1 0 3]; % [0 1;1 0] ambiguous
I(8,:)= [2 3 4 4]; % [1 1;1 0] 
I(9,:)= [3 2 4 4]; % [0 0;0 1] 
I(10,:)=[1 0 3 2]; % [1 0;0 1] ambiguous
I(11,:)=[3 1 4 4]; % [0 1;0 1] 
I(12,:)=[3 0 4 4]; % [1 1;0 1] 
I(13,:)=[0 2 4 4]; % [0 0;1 1] 
I(14,:)=[1 2 4 4]; % [1 0;1 1] 
I(15,:)=[0 1 4 4]; % [0 1;1 1] 
I(16,:)=[4 4 4 4]; % [1 1;1 1] 

% The base-edges by vertex positions 
E=[0 0 1 0; 0 0 0 1; 0 1 1 1; 1 0 1 1; 4 4 4 4];

% base-Polygons by vertexpostions 
IE=E(I(:)+1,:);
IE=[IE(1:16,:) IE(17:32,:); IE(33:48,:) IE(49:64,:)];

% Make a Binary image with pixels set to true above iso-treshold
B=Img>=isovalue;

% Get Elementary Cells
%  Cell
%  1 ** 2
%  *    *
%  *    *
%  4 ** 8
%
B0=B(1:end-1,1:end-1); B1=B(1:end-1,2:end); B2=B(2:end,1:end-1); B3=B(2:end,2:end);
V=B0+B1*2+B2*4+B3*8+1;
[x,y]=find((V>1)&(V<16));
v=V(x+(y-1)*size(V,1));

% Elementary cells to Edge coordinates defined by connected image grid-points
J=[IE(v,:);IE(v+16,:)];
r=J(:,1)==4;
J=J+[x y x y x y x y;x y x y x y x y];
J(r,:)=[];

% Vertices list defined by connected image grid-points
VP=[J(:,1:4);J(:,5:8)];

% Make a Face list
F=[(1:size(J,1))' (size(J,1)+1:2*size(J,1))'];

% Remove dubplicate vertices
[VP,a,Ind]=unique(VP,'rows'); F=Ind(F);

% Vertices described by image grid-points to real 
% linear Interpolated vertices
Vind1=VP(:,1)+(VP(:,2)-1)*size(Img,1);
Vind2=VP(:,3)+(VP(:,4)-1)*size(Img,1);
V1=abs(double(Img(Vind1))-double(isovalue));
V2=abs(double(Img(Vind2))-double(isovalue));
alpha=V2./(V1+V2);
Vx=VP(:,1).*alpha+VP(:,3).*(1-alpha);
Vy=VP(:,2).*alpha+VP(:,4).*(1-alpha);
V=[Vx Vy];





