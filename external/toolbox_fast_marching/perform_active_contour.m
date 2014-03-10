function D = perform_active_contour(D0, motion, options)

% perform_active_contour - perform active contour resolution
%
%   D = perform_active_contour(D0, motion, options);
%
%   Level set implementation of various active contour, all related to
%   motion by mean curvature.
%
%   The resolution is done by an explicit euler, the time sep is options.dt
%   and should be quite small for stability. The maximum time is
%   options.Tmax.
%
%   In the following, one denotes the curvature of the level sets of D as
%       Curv(D) = div(grad(D)/|grad(D)|)
%   and D' = d(D)/dt the time derivative.
%
%   motion can be:
%       'mean': D' = |grad(D)| * Curv(D)
%       'affine': D' = |grad(D)| * Curv(D)^(1/3)
%       'errosion': D' = |grad(D)| * max(Curv(D),0)^(1/3)
%       'snakes': D' = E * |grad(D)| * Curv(D) + <grad(E),grad(D)>
%       'chan-vese': D' = |grad(D)| * ( Curv(D) - lambda*(D-c1)^2 + lambda*(D-c2)^2 )
%
%   You can provide the parameters in options.E, options.c1, options.c2,
%   options.lambda. You can turn on the automatic update of c1/c2 using
%   options.update_c.
%
%   The distance function D is redistanced every options.redistance_freq
%   iterations.
%
%   To turn off the display, options.do_display=0. You can provide a
%   background image in options.M.
%
%   See also: display_segmentation.
%
%   Copyright (c) 2007 Gabriel eyre

dt              = getoptions(options, 'dt', 1000);
Tmax            = getoptions(options, 'Tmax', 1000);
do_display      = getoptions(options, 'display', 1);
display_freq    = getoptions(options, 'display_freq', 20);
redistance_freq = getoptions(options, 'redistance_freq', 30);
lambda          = getoptions(options, 'lambda', 0.8);
update_c        = getoptions(options, 'update_c', 1);
c1              = getoptions(options, 'c1', 0.1);
c2              = getoptions(options, 'c2', 0.8);
nb_svg          = getoptions(options, 'nb_svg', 0);
solver          = getoptions(options, 'solver', 'cg');
svg_path        = getoptions(options, 'svg_path', ''); % ['results/active-contour/' motion '/']);
if strcmp(motion,'chan-vese') || strcmp(motion, 'snake') || strcmp(motion, 'chan-vese-user')
    E           = getoptions(options, 'E', 0, 1);
else
    E = [];
end
M               = getoptions(options, 'M', E(:,:,1));
options.bound = 'per';

if not(isempty(svg_path)) && not(exist(svg_path))
    mkdir(svg_path);
end

nb_phase = 1;
if not(isempty(E))
    nb_phase = 2 - (size(E,3)==1);
end
    
if strcmp(motion, 'snake')
    % pre compute gradient of the energy
    dE = divgrad(E,options);
end

D = D0;
n = size(D0,1);
epsilon = 1e-3;
niter = round(Tmax/dt);

if nb_svg>0
    % distribute more svg at the begining of the iterations
    svg_list = round( 1 + linspace(0,1,nb_svg+2).^3*niter ); svg_list(1) = [];    
    % delta_svg = niter/nb_svg;
    if not(exist(svg_path))
        mkdir(svg_path);
    end
else
    delta_svg = Inf;
end
next_svg = 0;
nb = 0;


for i=1:niter
    progressbar(i,niter);
    for k=1:nb_phase
        g0(:,:,:,k) = divgrad(D(:,:,k),options);
        d(:,:,k) = max(epsilon, sqrt(sum(g0(:,:,:,k).^2,3)) );
        g(:,:,:,k) = g0(:,:,:,k) ./ repmat( d(:,:,k), [1 1 2] );
    end

    if strcmp(solver, 'grad')
        %%%% gradient %%%%
        switch lower(motion)
            case 'mean'
                dD = d .* divgrad( g,options );
            case 'affine'
                gg = divgrad(g,options);
                dD = d .* sign(gg) .* abs(gg).^(1/3);
            case 'errosion'
                dD = d .* max(divgrad(g,options),0).^(1/3);
            case 'snake'
                dD = E .* d .* divgrad( g,options ) + sum(dE.*g0, 3);
            case 'chan-vese'
                % compute the inner / outer constant
                if update_c
                    c1 = mean( E(D>=0) );
                    c2 = mean( E(D<0) );
                end
                dD = d .* ( divgrad( g,options ) - lambda*(E-c1).^2 + lambda*(E-c2).^2 );
            case 'chan-vese-user'
                if nb_phase==1
                    % single phase
                    dD = d .* ( divgrad( g,options ) + E );
                else
                    for k=1:nb_phase
                        dG = divgrad( g(:,:,:,k),options );
                        dD(:,:,k) = d(:,:,k) .* ( dG - E(:,:,2*k-1).*(D(:,:,3-k)>0) - E(:,:,2*k).*(D(:,:,3-k)<=0)  );
                    end                    
                end
            otherwise
                error('Unknown motion');
        end
        D = D + dt*dD;
    else
        options.E = [];
        switch lower(motion)
            case 'mean'
                C = zeros(n);
            case 'snake'
                options.E = E;
                C = sum(dE.*g0, 3);
            case 'chan-vese'
                % compute the inner / outer constant
                if update_c
                    c1 = mean( E(D>=0) );
                    c2 = mean( E(D<0) );
                end
                C = - lambda*(E-c1).^2 + lambda*(E-c2).^2;
            case 'chan-vese-user'
                C = E;
            otherwise
                error('Unknown motion');
        end        
        %%%% conjugate gradient %%%%
        % right hand side
        y = D + dt*C;
        options.d = d;
        options.dt = dt;
        options.n = n; options.ncols = n^2;
        options.niter_max = 20;
        ptions.epsilon = 1e-8;
        options.x = D(:);
        [D,err,k] = perform_conjugate_gradient(@callback_active_contour,y(:),options);
        D = reshape(D,n,n);
    end
        
    if mod(i,redistance_freq)==0
        for k=1:nb_phase
            D(:,:,k) = perform_redistancing(D(:,:,k), options);
        end
    end
    if do_display && ( mod(i,display_freq)==1 || i>=next_svg )
        A = prod(D,3);
        if strcmp(motion, 'snake') || strcmp(motion, 'chan-vese') || strcmp(motion, 'chan-vese-user')
            A = M;
        end
        clf;
        display_segmentation(D,A);
        drawnow;
        
        if 0
        clf;
        hold on;
        if strcmp(motion, 'snake') || strcmp(motion, 'chan-vese') || strcmp(motion, 'chan-vese-user')
            imagesc(M); axis image; axis off; axis([1 n 1 n]);
        else
            imagesc(D); axis image; axis off; axis([1 n 1 n]);
        end
%        h = plot(c(1,:), c(2,:), 'r');
        [c,h] = contour(D,[0 0], 'r');
        set(h, 'LineWidth', 2);
        drawnow;
        hold off;
        colormap gray(256);
        end
        
        if i>=next_svg
            nb = nb+1;
            % save image
            saveas(gcf, [svg_path motion '-' num2string_fixeddigit(nb, 2) '.png'], 'png');
            next_svg = svg_list(nb); % next_svg+delta_svg;
        end
    end
end