% test for the heuristically driven 
%   Fast Marching algorithm, 
%   aka FM*
%
%   Copyright (c) 2004 Gabriel Peyré

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speed fonction W

save_image = 1;

n = 600;

rep = 'images/fmstar_error/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

name = 'bump';
name = 'cavern';
name = 'map';
name = 'mountain';
name = 'cavern';
name = 'constant';
name = 'gaussian';
name = 'road2';

[M,W] = load_potential_map(name, n);
    

% pick start/end/center points
[start_points,end_points] = pick_start_end_point(M);
disp('Computing heuristic.');
[H,S] = perform_fast_marching(W, end_points);
disp('Computing reference distance');
[D,S] = perform_fast_marching(W, start_points);


weight_list = 0:0.1:1.1;

err = [];
i = 0;
clf;
for w = weight_list
    i = i+1;
    disp( sprintf('Propagation with weigth %.2f.', w) );
    options.null = 0;
    [D1,S] = perform_fast_marching(W, start_points, options, w*H);
    err = [err; norme(D-D1) ];
    subplot(3,3,min(i,9));
    imagesc(D-D1);
end

% plot error progression
clf;
plot( log10(weight_list(2:end)), log10(err(2:end)));
axis tight;
xlabel('log10(\lambda)');
ylabel('log10(error)');