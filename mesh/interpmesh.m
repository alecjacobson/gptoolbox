function [ d ] = interpmesh( V,F,ind,val,t)
%INTERPMESH Summary of this function goes here
%   Detailed explanation goes here

%% Build weight functions (initialization)

n = size(V,1);
nind = size(ind,2);

M = massmatrix(V,F,'voronoi');
L = cotmatrix(V,F);
a = sum(doublearea(V,F))/2.0;

W = zeros(nind,n);
% Build weight functions (computation)
for i=1:nind
    W(i,:) = prebiharmonic(M,L,a,t,ind(i));
end

meshplot(V,F,'ScalarFieldV',W(i,:)')

%% Compute weights for interpolation
c = W(:,ind)\val';

%% Interpolate function

d = zeros(n,1);
for i=1:nind
   d = d + c(i) * W(i,:)';
end

