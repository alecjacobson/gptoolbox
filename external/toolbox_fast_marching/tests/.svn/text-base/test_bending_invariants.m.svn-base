% computation of bending invariants


path(path, 'toolbox/');


disp('Loading mesh.');
% [vertex,faces] = gen_base_mesh('ico',4);

path(path, '../toolbox_graph_data/');
path(path, '../toolbox_graph_data/off/');

name = 'fandisk';
name = 'nefertiti';
name = 'beetle';
name = 'david50kf';
name = 'gargoyle';
name = 'horse';
name = 'bunny';
name = 'elephant-50kv';
name = 'hand';
[vertex,faces] = read_mesh(name);


options.nlandmarks = 200;
vertex1 = compute_bending_invariant(vertex,faces,options);


rep = ['results/bening-invariants/'];
if not(exist(rep))
    mkdir(rep);
end

% display both the original mesh and the embedding
options.name = name;
clf;
plot_mesh(vertex, faces, options);
shading interp; camlight;
saveas(gcf, [rep name '-original.png'], 'png');
if strcmp(name, 'david50kf') || strcmp(name, 'hand')
    zoom(.9);
end

clf;
plot_mesh(vertex1, faces, options);
shading interp; camlight;
saveas(gcf, [rep name '-invariant.png'], 'png');
if strcmp(name, 'david50kf') || strcmp(name, 'hand')
    zoom(.9);
end

