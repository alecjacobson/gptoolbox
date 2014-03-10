function [C] = connectedComponent(TT,TAG)
% CONNECTEDCOMPONENT Compute the connected components of a subset of triangles
% on a tri surface. Two faces are considered connected iff they share an edge.
%
%
% [C] = connectedComponent(TT,TAG)
%
% Inputs:
%   TT  face-face topology, Fx3 (you can use the function tt(F) to compute it)
%   TAG  subset of faces where the connected components are computed, Fx1
% Outputs:
%   C  cell array of connected components. Every element of the cell array
%      contains the indices of the faces in the corresponding connected
%      component
%
% Use size(C) to know the number of connected components detected
%

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
            for j=1:3
                neigh = TT(current,j);
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
