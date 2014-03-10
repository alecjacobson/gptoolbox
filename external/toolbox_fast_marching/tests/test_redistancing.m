% test for the re-distancing of a distance function
n = 256;

name = 'circle';
name = 'rectangle';
name = 'circlerect1';
name = 'circlerect2';

D0 = compute_levelset_shape(name, n);


% the modified distance function, should have the same 0 level set
D = (D0.^3) .* (X+n/3);

options.use_interpolation = 1;
D1 = perform_redistancing(D, options);
options.use_interpolation = 0;
D2 = perform_redistancing(D, options);

% original
c0 = perform_contour_extraction(D0, 0);
% interpolation
c1 = perform_contour_extraction(D1, 0);
% no interpolation
c2 = perform_contour_extraction(D2, 0);

clf;
hold on;
plot(c0(1,:),c0(2,:), 'r');
plot(c1(1,:),c1(2,:), 'g');
plot(c2(1,:),c2(2,:), 'b');
legend('Original', 'Interpolation', 'No interpolation');
axis([0 1 0 1]);
hold off;