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
