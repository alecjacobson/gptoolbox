function [S,SE] = sample_edges(V,E,samples_per_edge)
  % SAMPLE_EDGES Compute samples_per_edge extra points along each edge in E
  % defined over vertices of V.
  %
  % S = sample_edges(V,E,samples_per_edge)
  %
  % Inputs:
  %   V  vertices over which edges are defined, # vertices by dim
  %   E  edge list, # edges by 2
  %   samples_per_edge  number of extra samples to be computed along edge not
  %     including start and end points
  % Outputs:
  %   S  sampled vertices, size less than # edges * (2+samples_per_edge) by dim,
  %   always begins with V so that E is also defined over S
  %
  %

  dim = size(V,2);

  % trivial case
  if(isempty(E))
    S = [];
    SE = [];
  elseif(samples_per_edge < 0)
    S = V;
    SE = [];
  else
    % fraction parameter
    t = linspace(0.0,1.0,samples_per_edge+2);
    % get rid of start and end points
    t = t(2:(end-1));
    % repeat for each coordinate
    t = reshape(repmat(t,dim,1),1,size(t,2)*dim);
    % repeat for each edge
    t = repmat(t,size(E,1),1);
    % repeat start coords

    sp = repmat(V(E(:,1),:),1,samples_per_edge);
    % repeat end coords
    ep = repmat(V(E(:,2),:),1,samples_per_edge);
    % lerp from start point to end point for each coordinate
    S = sp.*(1-t) + ep.*t;
    % reshape to list coordinates
    S = S';
    S = reshape(S,dim,prod(size(S))/dim)';
    S = [V;S];
    % Determine edges between samples
    E1 = repmat((1:(samples_per_edge-1))',size(E,1),1);
    E2 = repmat((2:(samples_per_edge))',size(E,1),1);
    off = size(V,1) + ...
      reshape( ...
        repmat( ...
          samples_per_edge*((1:size(E,1))-1), ...
          samples_per_edge-1, ...
          1), ...
        (samples_per_edge-1)*size(E,1),...
        1);
    SE = repmat(off,1,2) + [E1 E2];
    % end point edges
    SE = [ ...
      SE; ...
      E(:,1) size(V,1) + (1+samples_per_edge*((1:size(E,1))-1))'; ...
      E(:,2) size(V,1) + (samples_per_edge*((1:size(E,1))))'];
  end
end
