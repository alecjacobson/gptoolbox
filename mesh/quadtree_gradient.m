function [G, stored_at] = quadtree_gradient(C,W,CH,D,A)
% QUADTREE_GRADIENT
% Builds a finite difference gradient on a quadtree following a centered 
% finite difference scheme, with the adjacency as suggested by 
% Bickel et al. "Adaptative Simulation of Electrical
% Discharges". This code is *purposefully* not optimized beyond
% asymptotics for simplicity in understanding its functionality and
% translating it to other programming languages beyond prototyping.
%
% G = quadtree_gradient(C,W,CH,D,A)
% [G,stored_at] = quadtree_gradient(C,W,CH,D,A)
%
% Inputs:
%   C #nodes by 3 matrix of cell centers
%   W #nodes vector of cell widths (**not** half widths)
%   CH #nodes by 4 matrix of child indeces (-1 if leaf node)
%   D #nodes vector of tree depths
%   A #nodes by #nodes sparse adjacency matrix, where a value of a in the
%       (i,j) entry means that node j is to the a-th direction of i
%       (a=1: left;  a=2: right;  a=3: bottom;  a=4: top).
%
% Outputs:
%   G #2*num_children by #num_children sparse gradient matrix (first
%       num_children rows are x derivatives, last are y derivatives)
%   stored_at #num_children by 3 matrix of child cell centers, where the
%       values of L are stored
%
%
% Example:
% %th = 2*pi*rand(200,1);
% P = 0.5*[cos(th),sin(th)];
% P = [P;[-1,-1];[1,1]];
% 
% % P = rand(500,2);
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MaxDepth',8,'Graded',true,'MinDepth',4);
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH);
% 
% [G, stored_at] = quadtree_gradient(C,W,CH,D,A);
% Gx = G(1:size(stored_at,1),:);
% Gy = G(size(stored_at,1) + (1:size(stored_at,1)),:);
% 
% % Sample function f(x,y)=x
% fun = stored_at(:,1);
% % This should be one almost everywhere
% disp(Gx*fun)
% hold off
% clf
% tsurf(Q,V,falpha(0,1))
% hold on
% sct(stored_at,30,'filled','CData',Gx*fun);
% set(gcf,'Color','w')
% grid off
% axis equal
% colorbar;
% % and it is!
% % This should be zero almost everywhere
% disp(Gy*fun)
% hold off
% clf
% tsurf(Q,V,falpha(0,1))
% hold on
% sct(stored_at,30,'filled','CData',Gy*fun);
% set(gcf,'Color','w')
% grid off
% axis equal
% colorbar;
% % and it... is close?
%
% See also: quadtree_laplacian.m



% We will store Laplacian values at
% child cell indeces
children = find(CH(:,1)==-1);
% map from all cells to children
cell_to_children = -ones(size(W,1),1);
cell_to_children(children) = 1:length(children);

% Vectors for constructing the Laplacian
I = [];
J = [];
vals = [];

for i=1:size(children,1)
    new_I = [];
    new_J = [];
    new_vals = [];
    l = [0,0,0,0];
    new_dirs = [];
    child = children(i);
    d = D(child);
    num_dirs = 0;
    % Let's build d u(child)/dx^2 ~ u(child+W(child)*[1,0])/hr(hl+hr) -
    % 2u(child)/hlhr + u(child-W(child)*[1,0])/hr(hl+hr)
    % So, let's look for the value to the j direction. To do this, we seek the
    % lowest-depth neighbor to the j direction. As a reminder the octree
    % adjacency convention is i->j (1:left-2:right-3:bottom-4:top)
    for j=1:4
        j_neighbors = find(A(child,:)==j);
        if ~isempty(j_neighbors)
            depths_j_neighbors = D(j_neighbors);
            [max_depth_j, max_depth_j_neighbor] = max(depths_j_neighbors);
            max_depth_j_neighbor = j_neighbors(max_depth_j_neighbor);
            % There are two options:
            % One: the leaf node to our j direction has lower or equal depth to
            % us
            if max_depth_j<=d
                l(j) = (W(child) + W(max_depth_j_neighbor))/2;
                % then it's easy, just add this node
                new_I = [new_I;i];
                % THIS HAS TO BE A CHILD !
                assert(cell_to_children(max_depth_j_neighbor)>0);
                new_J = [new_J;cell_to_children(max_depth_j_neighbor)];
                new_vals = [new_vals;1]; % Todo fix this
                new_dirs = [new_dirs;j];
            else
                % In this case, assuming the grid is graded, there should
                % be two j-neighbors at depth d+1
                nn = j_neighbors(D(j_neighbors)==(d+1));
                assert(length(nn)==2,"Are you sure you are inputting a graded quadtree?")
                assert(all(CH(nn,1)==[-1;-1]))
                % Then we simply average both
                l(j) = (W(child) + W(nn(1)))/2;
                new_I = [new_I;i;i];
                new_J = [new_J;cell_to_children(nn(1));cell_to_children(nn(2))];
                new_vals = [new_vals;0.5;0.5];
                new_dirs = [new_dirs;j;j];
            end
            num_dirs = num_dirs + 1;
        end
    end
    % This is a cheeky way to identify corners and make the stencil
    % backwards-forwards instead of centered in these cases
    for j=1:4
        if l(j)==0
            new_I = [new_I;i];
            new_J = [new_J;i];
            new_vals = [new_vals; 1.0];
            new_dirs = [new_dirs;j];
        end
    end
    % At this point, we have to divide by the edge-lengths and add sign
    new_vals(new_dirs==1) = -new_vals(new_dirs==1)/(l(1)+l(2));
    new_vals(new_dirs==2) = new_vals(new_dirs==2)/(l(1)+l(2));
    new_vals(new_dirs==3) = -new_vals(new_dirs==3)/(l(3)+l(4));
    new_vals(new_dirs==4) = new_vals(new_dirs==4)/(l(3)+l(4));
    % These are the y derivatives so they go in the lower matrix block
    new_I(new_dirs==3) = new_I(new_dirs==3) + size(children,1);
    new_I(new_dirs==4) = new_I(new_dirs==4) + size(children,1);

    % And add them to the big sparse Laplacian construction vectors
    I = [I;new_I];
    J = [J;new_J];
    vals = [vals;new_vals];
end

% THE LAPLACIAN IS NEGATIVE SEMI DEFINITE!
G = sparse(I,J,vals,2*length(children),length(children));
stored_at = C(children,:);

end