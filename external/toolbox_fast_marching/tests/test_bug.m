% addpath(genpath('fm'));

m = 10;
S = zeros(m,m,m); S(2:4,2:4,2:4)=1;

options.null = 0;

if 0
options.constraint_map = S;
options.constraint_map(S==0) = -Inf;
options.constraint_map(S~=0) = +Inf;
end


% gaussian weight (path will avoid center of the cube)
x = -1:2/(m-1):1;
[X,Y,Z] = meshgrid(x,x,x);
sigma = 0.4;
W = 1./(1 + exp( -(X.^2+Y.^2+Z.^2)/sigma^2 ) );

% options.nb_iter_max = Inf;

W = rescale(S,1e-5,1);
W = rand(m,m,m);
W = rescale(W,.5,1);


tic
[D,RS] = perform_fast_marching(W, [2; 2; 2], options);
fprintf(1,'m = %d, time = %f\n', m, toc);


return;

for m=[10 20 50 100]

    S = zeros(m,m,m); S(2:4,2:4,2:4)=1;

    options.constraint_map = S;
    options.constraint_map(S==0) = -Inf;
    options.constraint_map(S~=0) = +Inf;

    W = rescale(S,1e-5,1);

    tic
    [D,RS] = perform_fast_marching(W, [2; 2; 2], options);
    fprintf(1,'m = %d, time = %f\n', m, toc);
end


