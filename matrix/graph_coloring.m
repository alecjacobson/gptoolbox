function [C,nc] = graph_coloring(A,nc)
  % GRAPH_COLORING Color a graph given an adjacency matrix. This heuristic is
  % designed to efficiently find a conservative coloring (i.e., using too many
  % colors). It does *not* attempt to find the minimal coloring. It is best used
  % when the provided nc is a few more than the optimal coloring. For example,
  % planar graphs can be colored with 4 colors, (slow) heuristics for finding the
  % optimal coloring will often find 6 colors, but this function will very
  % quickly find a 7-coloring (nc=7). On such a graph, if you passed nc=4, this
  % function will likely waste time attemping to find a 4-coloring, 5-coloring,
  % 6-coloring, and then very quickly find a 7-coloring. It is designed to give
  % up searching for parsimonious colorings and increase the number of colors
  % (throwing a warning).
  %
  % [C,nc] = graph_coloring(A,nc)
  %
  % Inputs:
  %   A  #V by #V adjacency matrix
  %   nc  desired number of colors (for planar/triangle meshes, 7 is a fast
  %     choice; for tetrahedral meshes, 13 appears to often be fast)
  % Outputs:
  %   C  #V list of color ids into (1:nc)
  %   nc  effective number of colors >= input nc
  %
  % Examples:
  %   C = graph_coloring(facet_adjacency_matrix(F),9);
  %   tsurf(F,V,'CData',C);
  %   colormap(cbrewer('Set1',max(C)));

  n = size(A,1);
  [FI,FJ] = find(triu(A,1));
  % always marking the higher valence vertex as bad definitely makes things
  % worse 50x.
  %
  % if nc is much greater than optimal, always marking the lower valence vertex
  % doesn't much effect. But for nc close to optimal, this seems to make a very
  % big difference (300x)
  val = sum(A,2);
  swap = val(FI)<val(FJ);
  EI = FI.*swap + FJ.*~swap;
  EJ = FJ.*swap + FI.*~swap;


  onc = nc;
  % worst case just increase number of colors
  max_outer = 10;
  for outer = 1:max_outer
    % Might as well retry a few times for this nc
    for retry = 1:10
      C = ceil(nc*rand(n,1));
      for iter = 1:2*ceil(sqrt(size(A,1)))
        bad = unique(EI(C(EI)==C(EJ)));
        if numel(bad) == 0
          break;
        end
        C(bad) = nan;
        C(bad) = ceil(nc*rand(numel(bad),1));
      end
      if numel(bad) == 0
        break;
      end
    end
    if numel(bad) == 0
      break;
    end
    nc = nc+1;
    if nc>max_outer
      error('Failed to converge.\n');
    end
  end
  if nc ~= onc
    warning('Had to increasing number of colors')
  end

  %V2C = sparse(1:n,C,1,n,nc);
  % find all edges with the same color
  %A*V2C;
end
