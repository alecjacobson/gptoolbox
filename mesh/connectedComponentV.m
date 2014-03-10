function [C] = connectedComponentV(VV,TAG)
% CONNECTEDCOMPONENTV Compute the connected components of a subset of vertices
% on a tri surface. Two vertices are considered connected iff an edge connects
% them.
%
% [C] = connectedComponent(VV,TAG)
%
% Inputs:
%   VV  vertex-vertex topology, should be a sparse matrix as returned by the
%       function vv(V,F)
%   TAG  subset of vertices where the connected components are computed, Vx1
% Output:
%   C  cell array of connected components. Every element of the cell array
%      contains the indices of the vertices in the corresponding connected
%      component Use size(C) to know the number of connected components
%      detected

TOVISIT = 1;
VISITED = 2;

TAG(TAG ~= 0) = TOVISIT; % mark as to visit
C = cell(1,1);
Crow = 1;

I = find(TAG);

for i = I'
    if (TAG(i) == 1) % must be visited
        CElem = [i];
        toVisit = [i];
        TAG(i) = VISITED;
        
        while size(toVisit) ~= 0
            current = toVisit(1);
            toVisit = toVisit(2:end);
            ns = find(VV(:,current));
            for j=1:size(ns)
                neigh = ns(j);
                if (neigh == -1)
                    continue %% border
                end
                
                if (TAG(neigh) == TOVISIT)
                    TAG(neigh) = VISITED;
                    CElem = [CElem neigh];
                    toVisit = [toVisit neigh];
                end
            end
        end
        
        C{Crow} = CElem;
        Crow = Crow + 1;
    end
end

end
