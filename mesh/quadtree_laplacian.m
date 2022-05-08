function [L, stored_at] = quadtree_laplacian(C,W,CH,D,A)
% QUADTREE_LAPLACIAN
% Builds a finite difference Laplacian on a quadtree following the scheme
% suggested by Bickel et al. "Adaptative Simulation of Electrical
% Discharges". This code is *purposefully* not optimized beyond
% asymptotics for simplicity in understanding its functionality and
% translating it to other programming languages beyond prototyping.
%
% L = octree_laplacian(C,W,CH,D,A)
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
%   L #num_children by #num_children sparse Laplacian matrix
%   stored_at #num_children by 3 matrix of child cell centers, where the
%       values of L are stored
%
%
% Example:
%
% % Build an octree
% P = 0.5*[cos(th),sin(th)];
% P = [P;[-1,-1];[1,1]];
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MaxDepth',8,'Graded',true);
% % This is for plotting
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH);
% % Call function to construct Laplacian
% [L, stored_at] = quadtree_laplacian(C,W,CH,D,A);
% % Dummy Laplacian function
% gt_fun = stored_at(:,1).^2.0;
% laplacian_fun = 2 + 0.*stored_at(:,1);
% % Find boundary
% Vmin = min(stored_at,[],1);
% Vmax = max(stored_at,[],1);
% is_boundary = (stored_at(:,1)<=(Vmin(1)+0.2)) + (stored_at(:,2)<=(Vmin(2)+0.2)) + ...
%     (stored_at(:,1)>=(Vmax(1)-0.2)) + (stored_at(:,2)>=(Vmax(2)-0.2));
% bb = find(is_boundary);
% bc = gt_fun(bb);
% % Solve as energy minimization
% u = min_quad_with_fixed(0.5*L,-laplacian_fun,bb,bc);
% % Plot solution
% tsurf(Q,V,falpha(1,1),'FaceVertexCData',u)
% hold on
% sct(stored_at(bb,:)) % Visualize boundary conditions
% set(gcf,'Color','w'); grid off; axis equal; colorbar; caxis([0 1]);
%
%
% See also: initialize_quadtree.m





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
    l = [1,1,1,1];
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
                new_vals = [new_vals;-1]; % Todo fix this
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
                new_vals = [new_vals;-0.5;-0.5];
                new_dirs = [new_dirs;j;j];
            end
            num_dirs = num_dirs + 1;
        end
    end
    new_I = [new_I;i];
    new_J = [new_J;i];
    new_vals = [new_vals; 1.0];
    new_dirs = [new_dirs;5]; % just a hack
    % At this point, we have to divide by the edge-lengths
    new_vals(new_dirs==1) = new_vals(new_dirs==1)/(l(1)*(l(1)+l(2)));
    new_vals(new_dirs==2) = new_vals(new_dirs==2)/(l(2)*(l(1)+l(2)));
    new_vals(new_dirs==3) = new_vals(new_dirs==3)/(l(3)*(l(3)+l(4)));
    new_vals(new_dirs==4) = new_vals(new_dirs==4)/(l(4)*(l(3)+l(4)));
    new_vals(new_dirs==5) = 1/(l(1)*l(2)) + 1/(l(3)*l(4));
    
    % And add them to the big sparse Laplacian construction vectors
    I = [I;new_I];
    J = [J;new_J];
    vals = [vals;new_vals];
end

% THE LAPLACIAN IS NEGATIVE SEMI DEFINITE!
L = -2*sparse(I,J,vals,length(children),length(children));
stored_at = C(children,:);

end