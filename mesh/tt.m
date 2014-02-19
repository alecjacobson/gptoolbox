function [Fp, Fi] = tt(F)
  % TT Build a face adjacency data structure for a **manifold** triangle mesh.
  % From each face we can find where on its neighboring faces it's incident.
  %
  % [Fp, Fi] = tt(F)
  %
  % Input:
  %   F  list of face indices, #faces by 3
  % Output:
  %   Fp  #faces by 3, where Fp(i,j) tells the index of the neighboring triangle
  %     to the jth edge of the ith triangle in F. -1 if the jth edge of the ith
  %     triangle in F is a border edge. (jth edge refers to the edge opposite
  %     the jth vertex: so triangle in a triangle (a,b,c), the 1st edge is
  %     b-->c, the 2nd is c-->b and the 3rd is a-->b
  %   Fi  #faces by 3, where Fi(i,j) tells the position on the neighboring
  %     triangle to the jth edge of the ith triangle in F. -1 if the jth edge 
  %     if the ith triangle in F is a border edge. Uses the same indexing of
  %     positions of edges on a triangle as above.
  %
  % For example:
  %
  % F = [ 1 3 2;
  %       2 3 4];
  % [Fp, Fi] = tt(F);
  % Fp = 
  %     2    -1    -1
  %    -1    -1     1
  % Fi = 
  %     3    -1    -1
  %    -1    -1     1
  %
            
  if size(F,2) == 3

    % get list of edges (1st edges then 2nd edges then 3rd edges)
    E = [F(:,2) F(:,3); F(:,3) F(:,1); F(:,1) F(:,2)];

    % sparse adjacency matrix where (i,j) = k, k is the index of i-->j in edge
    % list if i-->j AND j-->i exist in edge list. This way is slightly faster
    % than building a proper adjacency list first.
    Ei = sparse(E(:,1),E(:,2),1:size(E,1));
    % adjacency list of edges as if edges represent undirected graph
    unadj = Ei>0;
    % "non-adjacency" list, (i,j) > 0 iff i-->j exists but j-->i does not exist
    nonadj = unadj-(unadj');
    adj = ((unadj + nonadj)==1).*Ei;

    % need to get mapping from sparse ordering to original order
    [ii jj si] = find(adj);
    % this is slightly faster than sorting the above by ii, which is the same
    [ii jj v] = find(adj');
    % build map from edges to their corresponding faces
    E2F = repmat(1:size(F,1),1,3)';
    % initialize adjacency map to -1
    Fp = -ones(size(E,1),1);
    Fp(si) = E2F(v);
    Fp = reshape(Fp,size(F));
    % build map from edges to edge positions
    I = reshape(repmat((1:3),size(F,1),1),1,3*size(F,1))';
    % use corresponding edges to find positions of correponding faces
    Fi = -ones(size(E,1),1);
    Fi(si) = I(v);
    Fi = reshape(Fi,size(F));
  else
    error('Not supported yet...');
  end

end

