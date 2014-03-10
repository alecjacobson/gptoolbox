% test for distance function compression

rep = 'images/landmark_fmstar/';
if exist(rep)~=7
    mkdir(rep);
end
rep_eps = [rep 'eps/'];
if exist(rep_eps)~=7
    mkdir(rep_eps);
end

n = 150;

save_image = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% potential creation
name = 'stephanodiscusniagarae';
name = 'map';
name = 'constant';
name = 'bump';
name = 'gaussian';
name = 'cavern';
name = 'mountain';
name = 'road2';

[M,W] = load_potential_map(name, n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% landmarks creation
landmark_init = 'farthest';
landmark_init = 'rand';
p = 30;
if strcmp(landmark_init, 'rand')
    landmarks = floor(rand(2,p)*n)+1;
else
    landmarks = perform_farthest_point_sampling( W, [], p );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute distance map
fprintf('Computing distance map ');
DL = zeros(n,n,p);
for k=1:p
    fprintf('.');
    DL(:,:,k) = perform_fast_marching(W, landmarks(:,k));
end
fprintf('\n');

% perform distance map estimation
DE = zeros(n,n,p);
for k=2:p
    DE(:,:,k) = compute_heuristic_landmark(DL(:,:,1:k-1),landmarks(:,k));
end

% details coefficients
DW = DL-DE;

% plot some details
marker_size = 20;
m = mmax(DW(:,:,2:end));
clf;
for k=1:9
    subplot(  3,3,k );
    hold on;
    M = abs(DW(:,:,k));
    M(1) = m;   % avoid scaling
    imagesc( M' );
    axis tight; axis off;
    plot(landmarks(1,k), landmarks(2,k), 'kx', 'MarkerSize', 10);
    plot(landmarks(1,1:k-1)', landmarks(2,1:k-1)', 'k.', 'MarkerSize', marker_size);
    hold off;
    colormap jet(256);
end

% some histograms
clf;
nb_bins = 30;
for k=1:16
    M = abs(DW(:,:,k));
    M1 = abs(DL(:,:,k));
    subplot(  4,4,k );
    hist([M(:),M1(:)],nb_bins);
    title(['point ' num2str(k)]);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform vector quantization


X = DL;





