function [ijk,D,uD,uV,uI,uJ] = lipschitz_octree(origin,h0,max_depth,udf,sdf)
  % [ijk,D,uD,uV,uI,uJ] = lipschitz_octree(origin,h0,max_depth,udf,sdf)
  %
  % Example:
  %   udf = @(P) sqrt(point_mesh_squared_distance(P,V,F));
  %   sdf = @(P) signed_distance(P,V,F,'SignedDistanceType','fwn');
  %   origin = min(V);
  %   h0 = max(max(V) - min(V));
  %   max_depth = 7;
  %   leaf_h = h0 / 2^max_depth;
  %   origin = 0.5*(max(V)+min(V)) - h0/2 - leaf_h;
  %   h0 = h0 + 2*leaf_h;
  %   [igk,D,uD,uV,uI,uJ] = lipschitz_octree(origin,h0,max_depth,udf,sdf);
  %   uJ = reshape(uJ,[],8);
  %   igl_uJ = uJ(:,[8     7     5     6     4     3     1     2]);
  %   [mV,mF] = marching_cubes_sparse(uD,uV,igl_uJ,0);

  function [unique_corner_positions,I,J,h] = unique_corners(origin,h0,ijk,depth)
    h = h0./(2.^depth);
    ijk_corners = ijk + permute(dec2bin(0:8-1,3)=='1',[3 2 1]);
    codes = sum(ijk_corners .* (2.^depth+1).^(0:2),2);
    [~,I,J] = unique(codes);
    all_ijk_corners = reshape(permute(ijk_corners,[1 3 2]),[],3);
    unique_ijk_corners = all_ijk_corners(I,:);
    unique_corner_positions = origin + h * unique_ijk_corners(:,[2 1 3]);
  end
  
  function [uU,unique_U,unique_corner_positions,I,J] = eval_on_corners(origin,h0,ijk,depth,f)
    %h = h0./(2.^depth);
    %corners = zeros(size(ijk,1),size(origin,2),8);
    %for i = 1:8
    %  ijk_this = ijk + (dec2bin(i-1,3)=='1');
    %  corners(:,:,i) = origin + h * ijk_this(:,[2 1 3]);
    %end
    %for c = 1:size(ijk,1)
    %  for i = 1:8
    %    fprintf("%d,%d : %g, %g, %g\n",c,i,corners(c,:,i));
    %  end
    %end
    %fprintf("\n");
    % U = reshape(udf(reshape(permute(corners,[1 3 2]),[],size(V,2))),size(ijk,1),8);
  
    [unique_corner_positions,I,J] = unique_corners(origin,h0,ijk,depth);
    unique_U = f(unique_corner_positions);
    uU = unique_U(J);
    uU = reshape(uU,size(ijk,1),8);
  end

  if nargin < 5
    sdf = [];
  end

  ijk = [0 0 0];
  count = 0;
  for depth = 0:max_depth
    h = h0./(2.^depth);
    
    % unused lipschitz constant
    lambda = 1;
    
    if depth == max_depth && ~isempty(sdf)
      empty = necessary;
      [D,uD,uV,uI,uJ] = eval_on_corners(origin,h0,ijk,depth,sdf);
  
      necessary = all(D>0,2) | all(D<0,2);
      ijk = ijk(~necessary,:);
    else
      % There's still some duplication between levels (1/8 of the nodes were
      % already computed on a previous level, and 1/8 of those were computed
      % previously etc.)
      % 
      % If the lipshitz bounds are deteremined with a faster unsigned distance
      % then maybe it doesn't matter so much since fwn signing is slower anyway.
      [U,uU] = eval_on_corners(origin,h0,ijk,depth,udf);
      count = count + size(uU,1);
  
      sufficient_1 = any(U > h*sqrt(3),2);
      sufficient_2 = all(U > h*sqrt(3)/2,2);
      necessary = true;
      empty = necessary & (sufficient_1 & sufficient_2);
      
      % only keep non-empty cells
      ijk_old = ijk;
      ijk = ijk(~empty,:);
      % split all non-empty cells into 8 children
      ijk_prev2 = ijk*2;
      ijk_next = [
        ijk_prev2 + [0 0 0]; ...
        ijk_prev2 + [0 0 1]; ...
        ijk_prev2 + [0 1 0]; ...
        ijk_prev2 + [0 1 1]; ...
        ijk_prev2 + [1 0 0]; ...
        ijk_prev2 + [1 0 1]; ...
        ijk_prev2 + [1 1 0]; ...
        ijk_prev2 + [1 1 1]; ...
      ];
    end
    
    
    %tsurf(F,V);
    %hold on;
    %if depth == max_depth
    %  draw_boxes(origin, h0, ijk, depth, falpha(0,1),'LineWidth',1);
    %else
    %  draw_boxes(origin, h0, ijk_old, depth, falpha(0,1),'LineWidth',1);
    %  draw_boxes(origin, h0, ijk_next, depth+1, falpha(0,1),'Edgecolor','r');
    %end
    %%sct(reshape(permute(corners,[1 3 2]),[],size(V,2)),'filled','SizeData',100,'CData',reshape(D,[],1));
    %hold off;
    %axis equal;
    %drawnow
  
    %fprintf('%-3.1f   = %10d / %-10d \n',...
    %  size(ijk,1)/(2^depth)^3, size(ijk,1),(2^depth)^3);
    % To right-align a number in printf use %10.1f
    % To left-align a number in printf use %-10.1f
    % 
    %fprintf('  %-3.1f = %10d / %-10d \n',...
    %  size(ijk,1)/(2^depth)^2, size(ijk,1),(2^depth)^2);
    
    if depth == max_depth
      break;
    end
  
    ijk = ijk_next;
  end
end
