% test for active contours
%
%   Copyright (c) 2007 Gabriel Peyre


options.bound = 'per';
n = 32;
d = ones(n); 
d = rand(n);
d = repmat(d,[1 1 2]);
D = zeros(n^2);
for i=1:n^2
    a = zeros(n); a(i)=1;
    a = divgrad(divgrad(a,options)./d,options);
    D(:,i) = a(:);
end

n = 150;

path(path, 'data/');
path(path, 'toolbox/');

motion = 'affine';
motion = 'errosion';
motion = 'chan-vese';
motion = 'mean';
motion = 'snake';

% load image for active contour
name = 'chan-vese';
name = 'disk';
name = 'brain';
M = rescale( sum( load_image(name, n), 3) );


options.solver = 'grad';
options.solver = 'cg';
    
%% load original shape
options.null = 0;
if strcmp(motion, 'snake')
    namec = 'square';
    options.width = 0.95*n;
elseif strcmp(motion, 'chan-vese')
    namec = 'small-disks';
else
    namec = 'circlerect1';
    namec = 'square';
    namec = 'circlerect2';
end
D0 = compute_levelset_shape(namec, n);

options.center = [0.15 0.15]*n;
options.radius = 0.1*n;
D1 = compute_levelset_shape('circle', n,options);


%% set up the parameters
options.Tmax = 1000;
options.redistance_freq = 150;
options.M = M;
options.dt = 0.1;
options.dt = 3;
options.display_freq = 5;
options.nb_svg = 100;
if strcmp(motion, 'snake')
    options.Tmax = 2000;
    options.redistance_freq = 30;
    options.dt = 1;
    % compute edge-based energy
    sigma = 4; % blurring size
    G = divgrad( perform_blurring(M,sigma) );
    G = sum( G.^2, 3);
    eta = 0.01;
    E = 1 ./ (eta + G);
    options.E = rescale(E, 0.7, 1);
elseif strcmp(motion, 'chan-vese')
    options.redistance_freq = 30;
    options.E = M;
    options.lambda = 0.8;
    options.update_c = 0;
    options.c1 = 0;
    options.c2 = 0.5;
    options.dt = 0.2;
    options.dt = 2;
end


D = perform_active_contour(D0, motion, options);