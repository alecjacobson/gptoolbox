function grad = compute_grad(M,options)

% compute_grad - compute the gradient of an image using central differences
%
% grad = compute_grad(M,options);
%
%   'options' is a structure:
%   - options.h is the sampling step size on both direction (default 1).
%   - options.h1 is the sampling step size on X direction (default 1).
%   - options.h2 is the sampling step size on Y direction (default 1).
%   - options.type is the kind of finite difference.
%       type==2 is fwd differences, ie.
%           y(i) = (x(i)-x(i-1))/h, with special
%           care at boundaries.
%       type==1 is forward differences bilinearly interpolated in the
%           middle of each pixel (be aware that you have a shift of 1/2 on X and Y for
%           the location of the gradient).
%       type==1 is backward differences bilinearly interpolated in the
%           middle of each pixel (be aware that you have a shift of -1/2 on X and Y for
%           the location of the gradient).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'h1')
    options.h1 = 1;
end
h1 = options.h1;
if ~isfield(options, 'h2')
    options.h2 = 1;
end
h2 = options.h2;

if isfield(options, 'h')
    h1 = options.h;
    h2 = options.h;
end

[n,p] = size(M);

if isfield(options, 'type')
    type = options.type;
else
    type = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new code, use faster 2D differences
if type==1
    % central differences on X
    D1 = [M(2:end,:);M(end,:)];
    D2 = [M(1,:);M(1:end-1,:)];
    grad(:,:,1) = (D1-D2)/(2*h1); 
    grad(1,:,1) = ( 4*M(2,:) - 3*M(1,:) - M(3,:) )/(2*h1);
    grad(end,:,1) = -( 4*M(end-1,:) - 3*M(end,:) - M(end-2,:) )/(2*h1);
    % central differences on Y
    D1 = [M(:,2:end),M(:,end)];
    D2 = [M(:,1),M(:,1:end-1)];
    grad(:,:,2) = (D1-D2)/(2*h2); 
    grad(:,1,2) = ( 4*M(:,2) - 3*M(:,1) - M(:,3) )/(2*h2);
    grad(:,end,2) = -( 4*M(:,end-1) - 3*M(:,end) - M(:,end-2) )/(2*h2);
elseif type==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(:,2:end),M(:,end)] )/2;
    % fwd differences on X
    D1 = [MM(2:end,:);MM(end,:)];
    D2 = MM;
    grad(:,:,1) = (D1-D2)/h1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on X
    MM = ( M + [M(2:end,:);M(end,:)] )/2;
    % fwd differences on Y
    D1 = [MM(:,2:end),MM(:,end)];
    D2 = MM;
    grad(:,:,2) = (D1-D2)/h2;
elseif type==3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(:,1),M(:,1:end-1)] )/2;
    % fwd differences on X
    D1 = MM;
    D2 = [MM(1,:);MM(1:end-1,:)];
    grad(:,:,1) = (D1-D2)/h;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % accumulate on Y
    MM = ( M + [M(1,:);M(1:end-1,:)] )/2;
    % fwd differences on Y
    D1 = MM;
    D2 = [MM(:,1),MM(:,1:end-1)];
    grad(:,:,2) = (D1-D2)/h;
else
    error('This kind of differences is not supported.');
end
 
return;

if type~=1
    % compute the difference in the center of each square
    h = zeros(n,p,2);
    
    for j=1:p-1
        h(:,j,1) = ( grad(:,j,1)+grad(:,j+1,1) )/2;
    end
    for i=1:n-1
        h(i,:,2) = ( grad(i,:,2)+grad(i+1,:,2) )/2;
    end
    
    if type==2          % fwd differences
        grad(1:n-1,1:p-1,:) = h(1:n-1,1:p-1,:);        
    elseif type==3      % bwd differences
        grad(2:n,2:p,:) = h(2:n,2:p,:);        
    end   
    
end