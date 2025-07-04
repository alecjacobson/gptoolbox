function [EF,EI,uE,EMAP] = edge_flaps(F, varargin)
  % [EF,EI,uE,EMAP] = Build edge flaps and topological edge relationships 
  % for a **manifold** triangle mesh.
  %
  % [EF, EI, uE, EMAP] = edge_flaps(F)
  % [EF, EI, uE, EMAP] = edge_flaps(F, 'CornerIndexing')
  % [EF, EI, uE, EMAP] = edge_flaps(F, 'CyclicIndexing')
  %
  % This function computes the "edge flaps" structure, describing for each unique
  % edge which triangles it belongs to, and its local edge index within each triangle.
  %
  % Inputs:
  %   F  #F by 3 list of face indices
  % Optional argument (mode):
  %   'CornerIndexing' (default) - Uses corner-based indexing convention: local edge
  %     indices are defined relative to the corner opposite each edge.
  %   'CyclicIndexing' - Uses edge-based cyclic indexing convention: local edge
  %     indices follow the order of edges (1→2, 2→3, 3→1) in each triangle.
  % Outputs:
  %   EF   #E by 2 matrix of "edge flaps". Each row corresponds to a unique edge.
  %        EF(e,1) is the index of the first adjacent face (triangle) that contains
  %        edge e in the direction (uE(e,1) → uE(e,2)). EF(e,2) is the index of the
  %        second adjacent face containing the same edge but in the opposite direction
  %        (uE(e,2) → uE(e,1)).
  %        If an edge is a boundary edge, then one of the entries in EF(e,:) will be 0.
  %   EI   #E by 2 matrix of local edge indices within the adjacent triangles
  %        corresponding to EF. Indices are 1-based and depend on the indexing mode.
  %   uE   #uE by 2 list of unique edges as vertex index pairs.
  %   EMAP #F by 3 matrix mapping each directed edge in F to its unique edge index in uE. 
  
  % Default mode
  mode = 'CornerIndexing';

  % Parse optional argument
  if ~isempty(varargin)
    if ischar(varargin{1}) || isstring(varargin{1})
      mode = varargin{1};
    else
      error('Optional argument must be a string: ''CornerIndexing'' or ''CyclicIndexing''.');
    end
  end

  E = [F(:,2:3);F(:,[3 1]);F(:,1:2)];
  sE = sort(E,2);
  [uE,~,EMAP] = unique(sE,'rows');
  EMAP = reshape(EMAP,size(F));
  I = (1:size(uE,1))';
  [B1,E1] = ismember(uE,E,'rows');
  [B2,E2] = ismember(uE,fliplr(E),'rows');
  J1 = mod(E1(B1)-1,size(F,1))+1;
  J2 = mod(E2(B2)-1,size(F,1))+1;
  C1 = floor((E1(B1)-1)/size(F,1))+1;
  C2 = floor((E2(B2)-1)/size(F,1))+1;
  I1 = I(B1);
  I2 = I(B2);
  K1 = repmat(1,numel(I1),1);
  K2 = repmat(2,numel(I2),1);
  EF = full(sparse([I1;I2],[K1;K2],[J1;J2],size(uE,1),2));
  EI = full(sparse([I1;I2],[K1;K2],[C1;C2],size(uE,1),2));

  % Apply cyclic reordering if requested
  if strcmpi(mode, 'CyclicIndexing')
    % Remap local edge indices
    i1 = EI == 1;
    i2 = EI == 2;
    i3 = EI == 3;
    EI(i3) = 1;
    EI(i1) = 2;
    EI(i2) = 3;
    % Reorder EMAP columns
    EMAP = [EMAP(:,3), EMAP(:,1), EMAP(:,2)];
  elseif ~strcmpi(mode, 'CornerIndexing')
    error('Unknown indexing mode. Use ''CornerIndexing'' or ''CyclicIndexing''.');
  end
end
