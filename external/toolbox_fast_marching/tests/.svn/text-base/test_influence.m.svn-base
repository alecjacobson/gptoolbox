% test for computation influence

n = 30;
W = ones(n);

x = -1:2/(n-1):1;
[Y,X] = meshgrid(x,x);
sigma = 0.4;
% W = 1./(1 + exp( -(X.^2+Y.^2)/sigma^2 ) );

start_points = [3;3];
end_points = [n-3;7];
nb_iter_max = Inf;

[Y,X] = meshgrid(0:1/(n-1):1, 0:1/(n-1):1);
H = sqrt( ( X-X(end_points(1), end_points(2)) ).^2 + ( Y-Y(end_points(1), end_points(2)) ).^2 );

[D,S,father] = perform_front_propagation_2d_slow(W,start_points,end_points,nb_iter_max);

% plot father relations
P = sub2ind(size(W), end_points(1), end_points(2));
M = zeros(n);
while ~isempty(P)
    % pop front
    p = P(1);
    [i,j] = ind2sub(size(W),p);
    P(1) = [];
    M(p) = 1;
    % add father
    f = father(i,j,:);
    if f(1)>0 && M(f(1))==0 && isempty(find(P==f(1)))
        P = [P,f(1)];
    end
    if f(2)>0 && M(f(2))==0 && isempty(find(P==f(2)))
        P = [P,f(2)];
    end  
end