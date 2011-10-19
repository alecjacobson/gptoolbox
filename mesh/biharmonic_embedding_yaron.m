function [V EG DD] = biharmonic_embedding_yaron(coord,triang)
%   Computes Biharmonic embedding V the Euclidean distance in that space 
%is the biharmonic distance.
%IMPORTANT: we use eigenvector approximation, make sure to use enough (eigno param).

%PARAMETER: this how many eigenvectors to use
eigno = 40; %%


% %read mesh from the off file
% fid=fopen(off_file);
% fgetl(fid);
% nos = fscanf(fid, '%d %d  %d', [3 1]);
% nopts = nos(1);
% notrg = nos(2);
% coord = fscanf(fid, '%g %g  %g', [3 nopts]);
% coord = coord';
% triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
% triang=triang';
% triang=triang(:,2:4)+1; %%we have added 1 because the vertex indices start from 0 in off format
% fclose(fid);

notrg=size(triang,1);
nopts=size(coord,1);

tic;
bir=[1 1 1];
trgarea=zeros(1,notrg);
for i=1:notrg
    trgx=[coord(triang(i,1),1) coord(triang(i,2),1) coord(triang(i,3),1)];
    trgy=[coord(triang(i,1),2) coord(triang(i,2),2) coord(triang(i,3),2)];
    trgz=[coord(triang(i,1),3) coord(triang(i,2),3) coord(triang(i,3),3)];
    aa=[trgx; trgy; bir];
    bb=[trgx; trgz; bir];
    cc=[trgy; trgz; bir];
    area=sqrt(det(aa)^2+det(bb)^2+det(cc)^2)/2;
    trgarea(i)=area;
end


%find the approximate voronoi area of each vertex
AM = zeros(nopts, 1);
for i=1:notrg
    AM(triang(i,1:3)) = AM(triang(i,1:3)) + trgarea(i)/3;
end

A = sparse(nopts, nopts);
%now construct the cotan laplacian
for i=1:notrg
    for ii=1:3
        for jj=(ii+1):3
            kk = 6 - ii - jj; % third vertex no
            v1 = triang(i,ii);
            v2 = triang(i,jj);
            v3 = triang(i,kk);
            e1 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v2,1) coord(v2,2) coord(v2,3)];
            e2 = [coord(v2,1) coord(v2,2) coord(v2,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            e3 = [coord(v1,1) coord(v1,2) coord(v1,3)] - [coord(v3,1) coord(v3,2) coord(v3,3)];
            cosa = e2* e3'/sqrt(sum(e2.^2)*sum(e3.^2));
            sina = sqrt(1 - cosa^2);
            cota = cosa/sina;
            w = 0.5*cota;
            A(v1, v1) = A(v1, v1) - w;
            A(v1, v2) = A(v1, v2) + w;
            A(v2, v2) = A(v2, v2) - w;
            A(v2, v1) = A(v2, v1) + w;
        end
    end
end

T = sparse([1:nopts], [1:nopts], (AM), nopts, nopts, nopts);

disp('calculating eigen-stuff...');
tic

%compute eigendecomposition until a tolerance is met
tol = 0.000001; %%tolerance
logtol = -log(tol);
devam = 1;
options.disp = 0;

if nargin < 3
    t = 1;
end
    
while devam
    tic;
    [V,D,flag] = eigs(A, T, eigno, 0, options);
    if flag ~= 0
        errrrrrrror=1
    end
    DD = -diag(D); %%eigens are made positive
    ratio = DD(eigno).^2;%/(t*DD(2));
    if ratio<logtol
        eigno = ceil(eigno*logtol/ratio);
        devam = 1;
    else
        devam = 0;
    end
end

% %normalize eigenvectors to unit length
% for i=1:eigno
%     V(:,i) = V(:,i)./norm(V(:,i));
% end
EG=V;


%get biharmonic embedding
for i=2:eigno
    V(:,i) = V(:,i)/(DD(i));%  exp(DD(i)/(2*t*DD(2)));
end

toc
disp('done!');
% 
% m = size(V,1);
% dists = sqrt(sum((ones(m,1)*V(src_ind,:) - V).^2,2));

% 
% %read the file with indices of points to find pairwise distances between
% if nargin < 4
%     sampleidx = (1:nopts)'; %%if no argument for index file, compute all pairs
% else
%     fid=fopen(idx_file);
%     sampleidx = fscanf(fid, '%d', [1 inf]);
%     sampleidx=sampleidx + 1; %%we have added 1 because the vertex indices start from 0
%     fclose(fid);
% end
% samplesize = length(sampleidx);
% 
% %compute distances
% V = V(sampleidx,:);
% savedists = pdist(V);
% dists = squareform(savedists);
% 
% %save the distances in to file
% fid = fopen(dist_file,'w');
% fprintf(fid,'%d\n',savedists);
% fclose(fid);
% 
% %save in vtk format for visual inspection 
% %use Paraview to view this format
% source_id = 1;
% ofid = fopen('control_diffusion.vtk','w');
% fprintf(ofid, '# vtk DataFile Version 3.0\n');
% fprintf(ofid,'vtk output\n');
% fprintf(ofid,'ASCII\n');
% fprintf(ofid,'DATASET POLYDATA\n');
% fprintf(ofid,'POINTS %d float\n', nopts);
% fprintf(ofid,'%g %g %g\n', coord');
% fprintf(ofid,'POLYGONS %d %d\n', notrg, 4*notrg)
% fprintf(ofid,'3 %d %d %d\n', triang'-1);
% fprintf(ofid,'\n');
% fprintf(ofid,'POINT_DATA %d\n', nopts);
% fprintf(ofid,'SCALARS distance_from float\n');
% fprintf(ofid,'LOOKUP_TABLE default\n');
% fprintf(ofid,'%g\n', dists(:,source_id));
% fclose(ofid);
% 
