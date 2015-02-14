function A = monotonicity_matrix(S,F,varargin)
  % MONOTONICITY_MATRIX  Computes a "monotonicity adjacency matrix" which
  % describes the directed acyclic "monotonicity graph" of "Topology-based
  % Smoothing of 2D Scalar Fields with C1-Continuity" by [Weinkauf, Gingold and
  % Sorkine 2010]
  %
  % A = monotonicity_matrix(S,F,'ParameterName',ParameterValue)
  %
  % Inputs:
  %   S  #V by 1 list of scalar values defined over (V,F)
  %   F  #F by 3 list of triangles indices
  %   Optional:
  %     'OneInOneOut' followed by bool controlling whether to limit graph to
  %       have only at least one graph edge coming and one out of each vertex,
  %       {false}
  %     'OneInOneOutPolicy' followed by:
  %        {'max'}  use max difference to add two edges per vertex
  %        'min'  use min difference to add two edges per vertex
  %        'strict' like max but then go through and greedily remove weakest
  %          edges that are superfluous (vertices on either side already have
  %          one in one out without this edge)
  %     'EnforceExtrema' followed by whether to explicitly enforce all extrema
  % Outputs:
  %   A  sparse #V by #V matrix where A(i,j) is zero unless j ∈ N(i) and S(j) >
  %     S(i). That is, A(i,j) reveals "parents" where a parent is connected by a
  %     directed edge that points down in terms of S. (A node may have multiple
  %     parents, but there will be no cycles). In that case A(i,j) =
  %     S(j)-S(i)
  %
  % See also: adjacency_matrix
  %

  % mandatory inputs
  n = size(S,1);
  enforce_extrema = true;

  % default options
  one_in_one_out = false;
  policy = 'max';
  % parse optional inputs
  ii = 1;
  while(ii <= numel(varargin))
    switch varargin{ii}
    case 'OneInOneOut'
      ii = ii + 1;
      assert(ii<=numel(varargin));
      one_in_one_out = varargin{ii};
    case 'OneInOneOutPolicy'
      ii = ii + 1;
      assert(ii<=numel(varargin));
      policy = varargin{ii};
    case 'EnforceExtrema'
      ii = ii + 1;
      assert(ii<=numel(varargin));
      enforce_extrema = varargin{ii};
    otherwise
      error(['Unknown input param: ''' varargin{ii} '''']);
    end
    ii = ii+1;
  end


  if size(F,2) > 2
    E = edges(F);
  else
    E = F;
  end
  % get both directions
  E = [E;fliplr(E)];

  A = cell(size(S,2),1);
  for ii = 1:size(S,2)
    A{ii} = sparse( ...
      E(:,1), ...
      E(:,2), ...
      (S(E(:,1),ii)<S(E(:,2),ii)).*(S(E(:,2),ii)-S(E(:,1),ii)), ...
      n,n);
    if one_in_one_out
      switch policy
      case 'max'
        [up_v, up] = max(A{ii},[],2);
        [down_v, down] = max(A{ii},[],1);
        A{ii} = sparse([1:n down(:)'],[up(:)' 1:n],[up_v(:)' down_v(:)'],n,n);
      case 'min'
        [up_v,up] = minnz(A{ii}');
        [down_v,down] = minnz(A{ii});
        A{ii} = sparse([1:n down(:)'],[up(:)' 1:n],[up_v(:)' down_v(:)'],n,n);
      case 'strict'
        warning('Dont remember if this actually works');
        % keep track of original minima and maxima since these will be
        % exceptions
        orig_min = find(sum(A{ii}~=0,1)==0);
        orig_max = find(sum(A{ii}~=0,2)==0)';
        down_count = sum(A{ii}~=0,1)';
        up_count = sum(A{ii}~=0,2);
        while(true)
          old_down_count = down_count;
          old_up_count = up_count;

          % remove minimum up edge
          [up_v,up] = minnz(A{ii}');
          % remove minimum down edge
          [down_v,down] = minnz(A{ii});

          % exit if nothing changes
          if all(old_down_count == down_count) && all(old_up_count == up_count)
            break;
          end
        end
        %% First use max to trim edges
        %[down_v, down] = max(A{ii},[],2);
        %[up_v, up] = max(A{ii},[],1);
        %A{ii} = sparse([1:n up(:)'],[down(:)' 1:n],[down_v(:)' up_v(:)'],n,n);
      otherwise
        error(['Unknown OneInOneOutPolicy: ''' policy '''']);
      end

      if enforce_extrema
        [~,xi] = local_max(F,S(:,ii));
        [~,ni] = local_min(F,S(:,ii));
        A{ii}(...
          sub2ind(size(A{ii}), ...
            [E(ismember(E(:,1),ni),1);E(ismember(E(:,1),xi),2)], ...
            [E(ismember(E(:,1),ni),2);E(ismember(E(:,1),xi),1)])) = ...
           S([E(ismember(E(:,1),ni),2);E(ismember(E(:,1),xi),1)]) - ...
           S([E(ismember(E(:,1),ni),1);E(ismember(E(:,1),xi),2)]);
      end

    end
  end

  % No matter what we need to enforce extrema


  if size(S,2) == 1;
    A = A{1};
  end
end
