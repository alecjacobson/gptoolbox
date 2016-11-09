function loop = outline_loop(F,varargin)
  % OUTLINE_LOOP
  %
  % loop = outline_loop(F)
  %  or 
  % loop = outline_loop(O)
  %
  % Inputs:
  %   F  #F by 3 list of triangles
  %   or 
  %   O  #O by 2 list of outline edges (needs to be oriented)
  %   Optional:
  %     'Unoriented' followed by true or false
  % Outputs:
  %   loop  #loop by 1 list of outline endpoints on single loop
  %
  %
  % Example:
  %   % list all loops
  %   O = outline(F);
  %   loops = {};
  %   while ~isempty(O)
  %     loops{end+1} = outline_loop(O);
  %     O = O(~all(ismember(O,loops{end}),2),:);
  %   end
  % 
  % Known bugs: Only works for single loop
  %

  if size(F,2) == 3 
    O = outline(F);
  else
    O = F;
  end

  unoriented = true;

  ii = 1;
  while ii < numel(varargin)
    switch varargin{ii}
    case 'Unoriented'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      unoriented = varargin{ii};
    otherwise
      error('Unsupported parameter: %s',varargin{ii});
    end
    ii = ii+1;
  end

  % Should only have single loop
  [~,IM] = remove_unreferenced((1:max(F(:)))',O);
  O = IM(O);
  RIM = sparse(IM,1,1:max(IM));
  loop = graphpred2path(sparse(O(2:end,1),1,O(2:end,2))',O(1,2));
  %if max(conncomp(adjacency_matrix(O))) > 1
  %    warning('Only first loop will be returned');
  %end

  %if unoriented
  %  A = adjacency_matrix(O(2:end,:));
  %  %[~, pred] = shortest_paths(A,O(1,2));
  %  %loop = path_from_pred(pred,O(1,1));
  %  [~, loop] = graphshortestpath(A,O(1,2));
  %  loop
  %else
  %  pred = full(sparse(1,O(:,2),O(:,1)));
  %  loop = path_from_pred(pred,1);
  %end

  loop = RIM(loop);
end
