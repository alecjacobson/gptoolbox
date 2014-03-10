clear all; close all; clc;
%--------------------------------------------------------------------------
library_directory = '../fmlib/';
addpath(library_directory);
data_directory   = '../Data/';
addpath(data_directory);
matlab_directory   = '../matlab/';
addpath(matlab_directory);
%%
%==========================================================================
n = 500;
% create 2D vector field
s = [n n 1];
U = randn(n,n,2);
sigma = (n/200)*40;
for it=1:10
    U = perform_vf_normalization( perform_blurring(U, sigma) );
end;
U = perform_vf_normalization( U );
% test for various degree of anisotropy
aniso_list = [.01 .05 .1 .2 .5 1];



%% test for progressive propagation
aniso = 0.1;
V = cat(3, -U(:,:,2), U(:,:,1)); % orthogonal vector
T = perform_tensor_recomp(U,V, ones(n),ones(n)*aniso );
plot_tensor_field(T);
source_points = pick_start_points(T(:,:,1,1));
h = [1;1];
%%
tic
[U, dUx, dUy, V, L] = fm2dAniso(h, T, source_points);
toc
%%
figure; imshow(U, []); colormap(jet);
hold on;
%plot_tensor_field(T);