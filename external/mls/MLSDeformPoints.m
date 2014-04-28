function MLSDeformPoints(v,p,type,outputName)
% MLSDEFORMPOINTS  Deform a set of point interactively

% Check params:
if nargin<2 || isempty(p)
    % Selecting the line case:
    if nargin>=2 && size(p,1)==4
        % Getting the segment extremes:
        fig=figure; plot(v(1,:),v(2,:),'k.');
        p=getpoints;
        close(fig);
        
        % Cropping the last odd element:
        if mod(size(p,2),2)==1
            p = p(:,1:end-1);
        end
        
        % Mixing the coords:
        p = [p(:,1:2:end);p(:,2:2:end)];
    else
        % Getting the pivot points:
        fig=figure; plot(v(1,:),v(2,:),'k.');
        p=getpoints;
        close(fig);
    end
end
if nargin<3 type='rigid'; end
if nargin<4 outputName='fv'; end

% Generating the mlsd:
if size(p,1)==4
    % Segments case:
    mlsd = MLSD2DlinesPrecompute(p,v,type);
else
    % Points case:
    mlsd = MLSD2DpointsPrecompute(p,v,type);
end

% Generate the figure:
fig = figure; hold on; axis equal;
plh = plot(mlsd.v(1,:),mlsd.v(2,:),'k.');
title(['Deformation type: ',mlsd.type]);

% Generating the set of moving points:
handles = [];
for i=1:size(mlsd.p,2)
    % Getting the api:
    if size(p,1)==4
        % Segments case:
        h = imline(gca,mlsd.p([1,3],i),mlsd.p([2,4],i));
    else
        % Points case:
        h = impoint(gca,mlsd.p(1,i),mlsd.p(2,i));
    end
    handles = [handles,h];
    api = iptgetapi(h);
    
    % Adding the callback:
    api.addNewPositionCallback(@Deform);
end

% Adding to the figure my data:
data.mlsd = mlsd;
data.q = mlsd.p;
data.handles = handles;
data.plh = plh;
data.outputName = outputName;
data.fv = v;
set(fig,'UserData',data);

% Adding the closing funciton:
set(fig,'CloseRequestFcn',@Closing);

% ------------------------ LOCAL FUNCTIONS ------------------------

% Saving the morphed image:
function Closing(varargin)

% Getting the current figure:
fig = varargin{1};

% Getting the data:
data = get(fig,'UserData');

% Saving the image:
assignin('base', data.outputName, data.fv);

% Closing:
closereq;

% -----------------------------------------------------------------

% The deformation function:
function Deform(np)

% Getting the figure:
fig = gcf;

% Getting the data:
data = get(fig,'UserData');

% Obtaining the actual positions:
if size(data.q,1)==4
    % Looking for the handler position:
    h = get(gco,'Parent');
    pos = find(data.handles==h);
    
    % Lines case:
    data.q(:,pos) = [np(1,:)';np(2,:)'];
else
    % Looking for the handler position:
    h = gco;
    pos = find(data.handles==h);

    % Points case:
    data.q(:,pos) = np(:);
end

% Deforming:
data.fv = MLSD2DTransform(data.mlsd,data.q);

% Plotting:
delete(data.plh);
data.plh = plot(data.fv(1,:),data.fv(2,:),'k.');

% Saving the data:
set(fig,'UserData',data);
