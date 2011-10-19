function [ IDS ] = samplemesh( V,F,n )
% [IDS] = samplemesh(V,F,n)
%
% Takes n sample over the surface
% There are no guarantees of the spacing of samples but they
% tend to be uniformly spaced
%
% Input:
%   V  vertex list
%   F  face list
%   n  number of samples
% Output:
%   IDS indices of the sampled vertices

IDS = zeros(n,1);
IDS(1) = 1;

D = zeros(size(V,1),n);
D(:,1) = sqrt(sum((repmat(V(1,:),size(V,1),1)-V(:,:)).^2,2));

for i=2:n
    S = sum(1.0./(D(:,1:(i-1)).^2),2);
    [~, j] = min(S);
    IDS(i) = j;
    D(:,i) = sqrt(sum((repmat(V(j,:),size(V,1),1)-V(:,:)).^2,2));
end

IDS = unique(IDS);
