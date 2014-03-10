function display_segmentation(B,M)

% display_segmentation - display a level set segmentation
%
%   display_segmentation(B,M);
%
%   B should be <0 inside the region of interest
%   M is a background image
%   
%   See also: perform_active_contour
%
%   Copyright (c) 2007 Gabriel Peyre

if min(B(:))>=0
    if max(B(:))>2
        % more that one phase
        B = cat(3, B>2, mod(B,2) );
    end
    for k=1:size(B,3)
        B(:,:,k) = rescale(B(:,:,k),-1,1);
    end
end

nb_phase = size(B,3);
n = size(B,1);

M = rescale(M);

if nb_phase==1
    M1 = M;
    M1(B<0) = M1(B<0)*.2;
    M = cat(3,M,M1,M1);
elseif nb_phase==2
    % 2 phases
    M1 = M; M2 = M;
    M1(B(:,:,1)<0) = M1(B(:,:,1)<0)*.2;
    M2(B(:,:,2)<0) = M2(B(:,:,2)<0)*.2;
    M = cat(3,M,M1,M2);
else
    error('Only 1/2 phase accepted');
end

hold on;
imagesc(M); axis image; axis off; axis([1 n 1 n]);
for k=1:nb_phase
    [c,h] = contour(B(:,:,k),[0 0], 'r');
    set(h, 'LineWidth', 2);
end
hold off;
colormap gray(256);
axis ij;