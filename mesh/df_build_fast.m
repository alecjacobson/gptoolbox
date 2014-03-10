function [DF] = df_build_fast(V,C)
% DF_BUILD_FAST Build a volumetric distance field and spatial index that can be
% used to speed up approximate NN queries in 3D This code uses the kdtree
% toolbox
%
% [DF] = df_build_fast(V,C)
%
% Input:
%   V  Coordinates of the points
%   C  Column vector that defines the number of cells for every
%      dimension
% Output:
%   DF  structure to be passed to df_query to execute queries

if ~exist('C','var')
    C=[10;10;10];
end

MIN = min(V,[],1);
MAX = max(V,[],1);

S = (MAX-MIN)./C';

D = zeros(C'+1);
N = zeros(C'+1);

%% Build kd-tree
tree = kdtree_build( V );


%% Queries
for x=1:size(D,1)
    progressbar(x,size(D,1));
    for y=1:size(D,2)
        for z=1:size(D,3)
            i = [x,y,z];
            p = MIN + (i-1).*S;
            
            N(x,y,z) = kdtree_nearest_neighbor(tree,p);
            D(x,y,z) = normrow(p-V(N(x,y,z),:));
        end
    end
end

kdtree_delete(tree);

% Prepare DF struct
DF.MIN = MIN;
DF.MAX = MAX;
DF.D = D;
DF.N = N;
DF.C = C;
DF.S = S;
end
