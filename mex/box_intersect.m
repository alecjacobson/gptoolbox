% BOX_INTERSECT Given sets of axis-aligne boxes determine which pairs overlap.
%
% I = box_intersect(A1,A2)
%
% Inputs:
%   A1  #A by dim list of minimum corners 
%   A2  #A by dim list of maximum corners 
% Outputs:
%   I  #I by 2 list of indices into A such that box I(i,1) of set A intersects
%     with box I(i,2) of set A.
% 
% I = box_intersect(A1,A2,B1,B2)
%
% Inputs:
%   A1  #A by dim list of minimum corners 
%   A2  #A by dim list of maximum corners 
%   B1  #B by dim list of minimum corners 
%   B2  #B by dim list of maximum corners 
% Outputs:
%   I  #I by 2 list of indices into A such that box I(i,1) of set A intersects
%     with box I(i,2) of set B.
%
% Example:
%   % Find all segment-segment intersections up to tolerance
%   [B1,B2] = box_each_element(V,E);
%   I = box_intersect(B1,B2);
%   % prune incident edges
%   keep = all(E(I(:,1),[1 1 2 2]) ~= E(I(:,2),[1 2 1 2]),2);
%   I = I(keep,:);
%   sqrD = segment_segment_squared_distance( ...
%     V(E(I(:,1),1),:), V(E(I(:,1),2),:), V(E(I(:,2),1),:), V(E(I(:,2),2),:));
%   I = I(sqrD<1e-7,:);
