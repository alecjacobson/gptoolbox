function [D1,Z] = compute_distance_landmark(start_points, DL, landmark, landmark_method);

% compute_distance_landmark - compute an heuristic using landmark points
%
%   [D1,Z] = compute_distance_landmark(start_points, DL, landmark,landmark_method);
%
%   DL(:,:,i) is the distance map to the ith landmark point.
%   'D1' is an approximation of the distance to 'start_point'.
%
%   Copyright (c) 2005 Gabriel Peyré


if nargin<4
    landmark_method = 'align';
end

if size(landmark,1)~=2
    landmark = landmark';
end

n = size(DL,1);
p = size(landmark,2);

if strcmp(landmark_method,'multiproxy0')
    landmark_method = 'align';
end

if length(landmark_method)>10 && strcmp(landmark_method(1:10), 'multiproxy')
    % order, aka number of ldm used
    order = str2num(landmark_method(11:end));
    order = min(order, size(DL,3));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESSING
    % distance between landmarks
    DD = zeros(p,p);
    for i=1:p
        for j=1:p
            DD(i,j) = min(DL(landmark(1,i),landmark(2,i),j), DL(landmark(1,j),landmark(2,j),i));
        end
    end
    % first get the k first distances
    I = zeros(n,n,order);
    DI = zeros(n,n,order);
    DL1 = DL; % to avoid destructing DL
    for i = 1:order
        [DI(:,:,i),I(:,:,i)] = min(DL1,[], 3);
        % put Inf values
        for x=1:n
            for y=1:n
                DL1(x,y,I(x,y,i)) = Inf;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PROCESSING
    D1 = zeros(n);
    Z = zeros(n);
    for i=1:order
        % ith landmark for start point
        Ii = I(start_points(1),start_points(2),i);
        for j=1:order
            % check I(start_point(1),start_point(2),i)<->I(:,:,j)
            Di = DD( Ii, I(:,:,j) ); Di = reshape(Di,n,n);
            Da = Di - DI(start_points(1),start_points(2),i) - DI(:,:,j);
            D1 = max( D1, Da );
            m = find(D1 == Da);
            Z(m) = i + (j-1)*order;
        end
    end
    return;
end

switch lower(landmark_method)

    case 'align'
        T = repmat(DL(start_points(1),start_points(2),:), [n,n,1] );
        [D1,Z] = max( abs(T-DL), [], 3 );
    
    case 'proxy'
        % find voronoi tesselation
        [V,I] = min(DL,[],3);
        start_points_I = I(start_points(1),start_points(2));
        % corresponding centers (ie quantization of the position)
        cx = reshape(landmark(1,I),n,n);
        cy = reshape(landmark(2,I),n,n);
        cxs = cx(start_points(1),start_points(2));
        cys = cy(start_points(1),start_points(2));
        % distance d(phi(start_point),phi(y))
        DD = DL( cxs,cys, I); 
        DD = reshape(DD,n,n);
        % symmetrize : distance d(phi(y),phi(start_point))
        DD1 = matrix_sampling_get(DL(:,:,start_points_I), [cx(:)';cy(:)']);
        DD1 = reshape(DD1,n,n);
%        DD = ( DD + DD1 )/2;   % pbm, this could break monotony
        DD = min(DD,DD1);
        % distance d(phi(x),x)
        Dx = V(start_points(1),start_points(2));
        % formula d(x,y) \approx d(phi(x),phi(y))-d(phi(x),x)-d(phi(y),y)
        D1 = max(DD-Dx-V,0);
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLD CODE
% compute the approximated distance
D1 = zeros(n);
Z = zeros(n);
for i=1:n
    for j=1:n
        [D1(i,j),Z(i,j)] = max( abs( DL(start_points(1),start_points(2),:)-DL(i,j,:) ) );
    end
end

