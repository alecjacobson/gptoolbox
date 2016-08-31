function [I,J] = stam_order(F,A)
  % STAM_ORDER For a regular patch of 13 faces and 12 vertices of a mesh (V,F)
  % around a "center" facet f determined by all its corners having valence 6,
  % determine the reordering of vertices, so that the reordered mesh
  % (V(I,:),J(F)) will fit Jos Stam's ordering in "Evaluation of Loop
  % Subdivision Surfaces". The facet F(f,:) is mapped to [4 7 8] in Stam's
  % patch.
  %
  % [I,J] = stam_order(F)
  % [I,J] = stam_order(F,A)
  %
  % Input:
  %   F  13 by 3 list of triangle indices into some V
  %   A  max(F) by max(F) adjacency matrix {[]}
  % Output:
  %   I  12 long list of indices such that V(I,:) are the reorderd vertices.
  %   J  12 long list of indices such that J(F) are the reordered faces,
  %     {J=full(sparse(I,1,1:12))}
  %
  % Example:
  %   % Canonical regular patch
  %   WV = [-1 0;0 -1;-1 1;0 0;1 -1;-1 2;0 1;1 0;2 -1;0 2;1 1;2 0];
  %   F = [ ...
  %     1 3 4;4 2 1;4 5 2; ...
  %     6 7 3;7 4 3;4 7 8;8 5 4;8 9 5; ...
  %     6 10 7;10 11 7;7 11 8;11 12 8;12 9 8];
  %   % Scramble patch vertex order and face order
  %   R = randperm(size(WV,1));
  %   G = R(F(randperm(end),:));
  %   U = full(sparse([R R],repmat([1 2],size(WV,1),1),WV));
  %   % Uncover order
  %   [I,J] = stam_order(G);
  %   % plot before and after
  %   subplot(1,2,1);
  %   tsurf(G,U,'VertexIndices',1); 
  %   set(gca,'YDir','reverse');
  %   subplot(1,2,2);
  %   tsurf(J(G),U(I,:),'VertexIndices',1); 
  %   set(gca,'YDir','reverse');
  %

  % There should be 13 face in a regular patch
  assert(size(F,1)==13,'F should have 13 faces to be a regular patch');

  if nargin<2 || isempty(A)
    A = adjacency_matrix(F);
  end
  assert(nnz(any(A))  == 12,'F should reference 12 vertices');
  % valences 
  C6 = sum(A,2) == 6;
  % Locate the center face:
  f = find(all(C6(F),2),1);
  if isempty(f) 
    error('F is not a regular patch');
  end
  % Stams 4,7,8 vertices
  I = zeros(12,1);
  I(4) = F(f,1);
  I(7) = F(f,2);
  I(8) = F(f,3);
  A(:,F(f,:)) = 0;
  I(5)  = find(A(I(4),:) & A(I(8),:));
  I(3)  = find(A(I(7),:) & A(I(4),:));
  I(11) = find(A(I(8),:) & A(I(7),:));
  A(:,I(I~=0)) = 0;
  I(12) = find(A(I(11),:) & A(I(8),:));
  I( 9) = find(A(I(12),:) & A(I(8),:));
  I( 2) = find(A(I( 5),:) & A(I(4),:));
  I( 1) = find(A(I( 2),:) & A(I(4),:));
  I( 6) = find(A(I( 3),:) & A(I(7),:));
  I(10) = find(A(I( 6),:) & A(I(7),:));
  J = full(sparse(I,1,1:12));
end
