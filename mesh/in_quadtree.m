function [i, others] = in_quadtree(point,C,W,CH)
% IN_QUADTREE
% Traverses a quadtree in logarithmic time to find the smallest cell
% containing a given point in 2D space
%
% [i, others] = in_quadtree(point,C,W,CH)
%
% Inputs:
%   point size-3 vector of point in the plane
%   C #nodes by 3 matrix of cell centers
%   W #nodes vector of cell widths (**not** half widths)
%   CH #nodes by 4 matrix of child indeces (-1 if leaf node)
%
% Outputs:
%   i integer index of smallest cell containint P into C,W,CH
%   others vector of integers to all other (non-leaf) cells containing P
%
% Example:
%
% th = 2*pi*rand(200,1);
% P = 0.5*[cos(th),sin(th)];
% P = [P;[-1,-1];[1,1]];
% 
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MaxDepth',8,'Graded',true);
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH);
% hold on
% t = tsurf(Q,V,'FaceColor','b','EdgeColor','k',falpha(0.2,1));
% cage = C(1,:) + 0.5*W(1)*[-1,-1;1,-1;1,1;-1,1;-1,-1];
% r = plot(cage(:,1),cage(:,2),'-r','LineWidth',3);
% axis equal
% 
% queries = 2*rand(100,2)-1;
% q = plot(0,0,'.k','MarkerSize',30);
% 
% for i=1:size(queries,1)
%     q.XData = queries(i,1);
%     q.YData = queries(i,2);
%     ind = in_quadtree(queries(i,:),C,W,CH);
%     cage = C(ind,:) + 0.5*W(ind)*[-1,-1;1,-1;1,1;-1,1;-1,-1];
%     r.XData = cage(:,1);
%     r.YData = cage(:,2);
%     figgif('queries.gif')
% end
%
%
%
others = [];
queue = 1;
i = -1; % by default it's nowhere

function b = is_in_quad(queries,center,width)
    max_corner = center + width*[0.5,0.5];
    min_corner = center - width*[0.5,0.5];
    b = ( (queries(:,1)>=min_corner(1)) & (queries(:,2)>=min_corner(2)) ...
        & (queries(:,1)<=max_corner(1)) & (queries(:,2)<=max_corner(2)) );
end

while ~isempty(queue)
    % Pop from queue
    q = queue(1);
    queue(1) = [];
    % Check if point is inside this cell
    if is_in_quad(point,C(q,:),W(q))
        % If inside this cell, then add to otthers
        others = [others;q];
        % Is it child?
        is_child = (CH(q,1)==-1);
        if is_child
            % If it is, then we're done
            i = q;
            break
        else
            % If not, add children to queue
            queue = [queue;CH(q,:)'];
        end
    end
end

end