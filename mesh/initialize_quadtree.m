function [C,W,CH,PAR,D,A] = initialize_quadtree(P,varargin)
% INITIALIZE QUADTREE
% Builds an adaptatively refined (optionally graded) quadtree for
% prototyping on adaptative grids. Keeps track of all parenthood and
% adjacency information so that traversals and differential quantities are
% easy to compute. This code is *purposefully* not optimized beyond
% asymptotics for simplicity in understanding its functionality and
% translating it to other programming languages beyond prototyping.
%
% 
% [C,W,CH,PAR,D,A] = initialize_quadtree(P)
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MinDepth',3,'MaxDepth',9)
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'Graded',true)
%
%
% Inputs:
%   P is a #P by 3 matrix of points. The output tree will be more subdivided in
%       regions with more points
%   Optional:
%       MinDepth integer minimum tree depth (depth one is a single box)
%       MaxDepth integer max tree depth (min edge length will be
%           bounding_box_length*2^(-MaxDepth))
%       Graded boolean whether to ensure that adjacent quads only differ by
%           one in depth or not (this is useful for numerical applications, 
%           not so much for others like position queries).
%
% Outputs:
%   C #nodes by 3 matrix of cell centers
%   W #nodes vector of cell widths (**not** half widths)
%   CH #nodes by 4 matrix of child indeces (-1 if leaf node)
%   PAR #nodes vector of immediate parent indeces (to traverse upwards)
%   D #nodes vector of tree depths
%   A #nodes by #nodes sparse adjacency matrix, where a value of a in the
%       (i,j) entry means that node j is to the a-th direction of i
%       (a=1: left;  a=2: right;  a=3: bottom;  a=4: top).
%
%
% Example:
% th = 2*pi*rand(200,1);
% P = 0.5*[cos(th),sin(th)];
% P = [P;[-1,-1];[1,1]];
% [C,W,CH,PAR,D,A] = initialize_quadtree(P,'MaxDepth',8,'Graded',true);
% [V,Q] = bad_quad_mesh_from_quadtree(C,W,CH);
% 
% clf
% hold off
% hold on
% tsurf(Q,V,'FaceColor','b','EdgeColor','k',falpha(0.2,1))
% axis equal
%
% See also: bad_quad_mesh_from_quadtree.m


% default values
min_depth = 1;
max_depth = 7;
graded = true;
% Map of parameter names to variable names
params_to_variables = containers.Map( ...
{'MaxDepth','MinDepth','Graded'}, ...
{'max_depth','min_depth','graded'});
v = 1;
while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
end


% Auxiliary function
function R = transpose_orientation(L)
    R = L';
    % left to right
    R(L'==1) = 2;
    R(L'==2) = 1;
    % top to bottom
    R(L'==3) = 4;
    R(L'==4) = 3;
end

% Check if point in box
function b = is_in_quad(queries,center,width)
    max_corner = center + width*[0.5,0.5];
    min_corner = center - width*[0.5,0.5];
    b = ( (queries(:,1)>=min_corner(1)) & (queries(:,2)>=min_corner(2)) ...
        & (queries(:,1)<=max_corner(1)) & (queries(:,2)<=max_corner(2)) );
end

% Subdivide
function [C2,W2,CH2,PAR2,D2,A2] = subdivide_quad(ind,C1,W1,CH1,PAR1,D1,A1,graded1)
    % Can't subdivide something that is not currently a leaf node
    assert(CH1(ind,1)==-1)
    % For simplicity:
    w = W1(ind);
    c = C1(ind,:);
    d = D1(ind);
    p = PAR1(ind);
    a_ind = A1(:,ind); % from the perspective of the neighbors
    num_quads = size(C1,1);
    % Easy part: Add four new cells
    % order: bottom-left, bottom-right, top-left, top-right
    C2 = [C1;...
        c + 0.25*w*[-1,-1];...
        c + 0.25*w*[1,-1];...
        c + 0.25*w*[-1,1];...
        c + 0.25*w*[1,1]];
    % New cells have half the length
    W2 = [W1;0.5*w;0.5*w;0.5*w;0.5*w];
    % New cells have one more depth
    D2 = [D1;d+1;d+1;d+1;d+1];
    % Keep track of child indeces
    CH2 = [CH1;repmat([-1,-1,-1,-1],4,1)];
    CH2(ind,:) = num_quads + [1,2,3,4];
    % and parent indeces
    PAR2 = [PAR1;ind;ind;ind;ind];
    % Now the hard part, which is the adjacency
    % (left-right-bottom-top)
    % Effectively we are concatenating [A , B;  "-B", C]
    % C is always the same 4 by 4 matrix
    square_mat = sparse([0,2,4,0;1,0,0,4;3,0,0,2;0,3,1,0]);
    % Now the B matrix is num_quads by 4
    rect_mat = sparse(num_quads,4);
    % We populate it by going over the neighbors of ind
    neighbors_ind = find(a_ind);
    where_neighbors = a_ind(neighbors_ind);
    % These should all be non-zero
    for i=1:length(neighbors_ind)
        neighbor_ind = neighbors_ind(i);
        neighbor_where = where_neighbors(i);
        neighbor_depth = D1(neighbor_ind);
        % This will help us populate rect_mat(:,neighbor_ind)
        % We'll build I, J and val.
        % There are two options here: if the neighbor has the same or
        % higher depth, it will gain two new neighbors. If not, it will
        % gain only one. Maybe there's a way to combine both options but
        % for now let's do it one at a time.
        
        % Let's start with the easy case: neighbor depth is low
        if neighbor_depth<=d
            % J will always be the same (neighbor_ind will gain two
            % neighbors)
            J = [neighbor_ind,neighbor_ind];
            % Orientation will be the same:
            vals = [neighbor_where, neighbor_where];
            % The tricky bit is *which* are the new neighbors
            % order: bottom-left, bottom-right, top-left, top-right
            if neighbor_where==1 % if ind_quad is to the left of neighbor_ind
                I = [2,4]; % right indeces
            elseif neighbor_where==2 % if ind_quad is to the right of neighbor_ind
                I = [1,3]; % left indeces
            elseif neighbor_where==3 % if ind_quad is to the bottom of neighbor_ind
                I = [3,4]; % top indeces
            elseif neighbor_where==4 % if ind_quad is to the top of neighbor_ind
                I = [1,2]; % bottom indeces
            end
        else
            % neighbor depth is high. We need to traverse the tree to find
            % out which of the four d + 1 depth children this neighbor
            % comes from originally. This, combined with the neighbor_where
            % information, will tell us which one new neighbor this cell
            % has. Should use a look-up table with 16 possibilities
            % J will always be the same (neighbor_ind will gain 1 neighbor)
            J = [neighbor_ind];
            % Orientation will be the same:
            vals = [neighbor_where];
            % 
            % First,  which of the four depth d + 1 children do we come
            % from
            n_ind = neighbor_ind;
            n_depth = neighbor_depth;
            n_par = PAR1(neighbor_ind);
            while n_depth>d+1 % it may be that we don't ever enter this loop
                n_ind = n_par;
                assert(n_depth==(D1(n_ind)+1))
                n_depth = D1(n_ind);
                n_par = PAR1(n_ind);
            end
            which_child = find(CH1(n_par,:)==n_ind);
            assert(d==D1(n_par))
            % Reminder: bottom-left, bottom-right, top-left, top-right
            if which_child==1 % it comes from the bottom left bit of the depth-d neighbor
                % Then there are two options, either ind_quad is to its
                % left or to its bottom
                if neighbor_where==1 % if ind_quad is to the left of neighbor_ind
                    I = 2; % then this is the bottom right of neighbor_ind
                elseif neighbor_where==3 % if ind_quad is to the bottom of neighbor_ind
                    I = 3; % then this is the top left
                end
            elseif which_child==2 % it comes from the BOTTOM RIGHT bit of the depth-d neighbor
                % Then there are two options, either ind_quad is to its
                % right or to its bottom
                if neighbor_where==2 % right
                    I = 1; % bottom left
                elseif neighbor_where==3 % bottom
                    I = 4; % top right
                end
            elseif which_child==3 % it comes from the TOP LEFT bit of the depth-d neighbor
                % two option: top or left
                if neighbor_where==1 % left
                    I = 4; % top right
                elseif neighbor_where==4 % top
                    I = 1; % bottom left
                end
            elseif which_child==4 % it comes from the TOP RIGHT bit
                % two options: top or right
                if neighbor_where==4 % top
                    I = 2; % bottom right
                elseif neighbor_where==2
                    I = 3;
                end
            end
            % ...phew. Again, that should be a lookup table
        end
        rect_mat = rect_mat + sparse(J,I,vals,num_quads,4);
    end
    A2 = [A1,rect_mat;transpose_orientation(rect_mat),square_mat];
    % Shouldn't matter but just for consistency let's use A2, ...2 info
    if graded1
        a2_ind = A2(:,ind); % from the perspective of the neighbors
        % Go over all neighbors of the original quad we subdivided
        neighbors_ind = find(a2_ind);
        for i=1:length(neighbors_ind)
            is_child_2 = ((CH2(neighbors_ind(i),1)==-1));
            neighbor_depth = D2(neighbors_ind(i));
            if (is_child_2 && neighbor_depth<d)
                [C2,W2,CH2,PAR2,D2,A2] = subdivide_quad(neighbors_ind(i),C2,W2,CH2,PAR2,D2,A2,true);
            end
        end
        
    end
end



% We start with a bounding box
Vmin = min(P,[],1);
Vmax = max(P,[],1);
C = (Vmin + Vmax)/2;
W = max(Vmax-Vmin);
CH = [-1,-1,-1,-1]; % for now it's leaf node
D = 1;
A = sparse(1,1);
PAR = -1; % supreme Neanderthal ancestral node


% Now, we will loop
quad_ind = 0;
while true
    quad_ind = quad_ind + 1;
    if quad_ind>size(C,1)
        break
    end
    is_child = (CH(quad_ind,1)==-1);
    % Does this quad contain any point? (Or is it below our min depth)
    if ((D(quad_ind)<min_depth || any(is_in_quad(P,C(quad_ind,:),W(quad_ind)))) ...
            && D(quad_ind)<max_depth && is_child)
        % If it does, subdivide it
        [C,W,CH,PAR,D,A] = subdivide_quad(quad_ind,C,W,CH,PAR,D,A,graded);
    end
end


end