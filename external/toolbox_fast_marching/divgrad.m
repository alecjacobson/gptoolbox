function G = divgrad(M,options)

% divgrad - compute either gradient or divergence.
%
%   G = divgrad(M);
%
%   if M is a 2D array, compute gradient, 
%   if M is a 3D array, compute divergence.
%   Use centered finite differences.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if size(M,3)==2
    G = mydiv(M,options);
else
    G = mygrad(M,options);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = mydiv(g,options)

bound = getoptions(options, 'bound', 'sym');

n = size(g,1);
d = zeros(n);

if strcmp(bound,'sym')
    d(2:end-1,:) = d(2:end-1,:) + ( g(3:end,:,1)-g(1:end-2,:,1) )/2;
    d(1,:) = d(1,:) + g(2,:,1)-g(1,:,1);
    d(end,:) = d(end,:) + g(end,:,1)-g(end-1,:,1);

    d(:,2:end-1) = d(:,2:end-1) + ( g(:,3:end,2)-g(:,1:end-2,2) )/2;
    d(:,1) = d(:,1) + g(:,2,2)-g(:,1,2);
    d(:,end) = d(:,end) + g(:,end,2)-g(:,end-1,2);
else
    sel1 = [2:n 1];
    sel2 = [n 1:n-1];
    d = g(sel1,:,1)-g(sel2,:,1) + g(:,sel1,2)-g(:,sel2,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function g = mygrad(M,options)

bound = getoptions(options, 'bound', 'sym');

n = size(M,1);
g = zeros(n,n,2);

if strcmp(bound,'sym')
    % on x
    g(2:end-1,:,1) = ( M(3:end,:)-M(1:end-2,:) )/2;
    g(1,:,1) = M(2,:)-M(1,:);
    g(end,:,1) = M(end,:)-M(end-1,:);
    % on y
    g(:,2:end-1,2) = ( M(:,3:end)-M(:,1:end-2,:) )/2;
    g(:,1,1) = M(2,:)-M(1,:);
    g(:,end,1) = M(:,end)-M(:,end-1);
else
    sel1 = [2:n 1];
    sel2 = [n 1:n-1];
    g = cat( 3, M(sel1,:)-M(sel2,:), M(:,sel1)-M(:,sel2) )/2;
end