function [VV,FF,J] = insert_vertices(V,F,FI,B)
  % [VV,FF,J] = insert_vertices(V,F,FI,B)
  %
  % Insert vertices into a mesh defined by vertices V and faces F
  % at the barycentric coordinates B in the faces FI.
  %
  % Inputs:
  %   V   #V by 3 list of vertex positions
  %   F   #F by 3 list of triangle indices into V
  %   FI  #P by 1 list of face indices into F
  %   B   #P by 3 barycentric coordinates in corresponding face in FI
  % Outputs:
  %   VV  #VV by 3 list of vertex positions after insertion
  %   FF  #FF by 3 list of triangle indices into VV
  %   J   #FF by 1 list of face indices into F

  assert(numel(FI) == size(B,1));
  assert(all(abs(sum(B,2)-1)<1e-10));

  J = (1:size(F,1))';
  K = zeros(numel(FI),1);

  NZ = sum(B==0,2);
  % Vertices
  IV = NZ==2;
  % Edges
  IE = NZ==1;
  % Faces
  IF = NZ==0;
  assert(all(IV+IE+IF == 1));

  % Vertices are no-ops
  FIV = FI(IV);
  BV = B(IV,:);
  [BVI,BVJ,BVV] = find(BV);
  assert(all(BVV)==1);
  K(IV(BVI)) = F(sub2ind(size(F),FI(BVI),BVJ));
  VV = V;
  FF = F;

  % inserting in a face only affects that face
  FIF = FI(IF);
  FI2F = sparse(1:numel(FIF),FIF,1,numel(FIF),size(F,1));
  BF = B(IF,:);

  A = FI2F*FI2F';
  %warning('rng');
  %rng(6);
  C = graph_coloring(A,max(sum(FI2F)));
  for c = reshape(unique(C),1,[])
    pop = full(C==c);
    FIFpop = FIF(pop);
    Vpop =  ...
     V(FF(FIFpop,1),:).*BF(pop,1) + ...
     V(FF(FIFpop,2),:).*BF(pop,2) + ...
     V(FF(FIFpop,3),:).*BF(pop,3);
    npop = sum(pop);
    % Split Face
    % Append new faces
    m = size(FF,1);
    S1 = m+0*npop+(1:npop)';
    S2 = m+1*npop+(1:npop)';
    S3 = m+2*npop+(1:npop)';

    A = sparse(size(F,1),3);
    A(FIFpop,:) = BF(pop,:);
    A = A(FI,:);
    BdA = B./A;
    G = zeros(size(B,1),3,3);
    for d = 1:3
      d1 = mod(d-1+0,3)+1;
      d2 = mod(d-1+1,3)+1;
      d3 = mod(d-1+2,3)+1;
      G(:,d1,d1) = BdA(:,d1);
      G(:,d2,d1) = B(:,d2)-BdA(:,d1).*A(:,d2);
      G(:,d3,d1) = B(:,d3)-BdA(:,d1).*A(:,d3);
    end
    [minG,D] = max(min(G,[],2),[],3);
    G1 = G(:,:,1);
    G2 = G(:,:,2);
    G3 = G(:,:,3);
    G = zeros(size(B));
    G(D==1,:) = G1(D==1,:);
    G(D==2,:) = G2(D==2,:);
    G(D==3,:) = G3(D==3,:);
    GI = zeros(size(B,1),1);
    npop
    size(D)
    GI(D==1,:) = S1(D==1);
    GI(D==2,:) = S2(D==2);
    GI(D==3,:) = S3(D==3);
    assert(all(minG>0));
    error



    AI = size(VV,1)+(1:npop)';
    FF = [FF
       AI FF(FIFpop,2) FF(FIFpop,3)
       FF(FIFpop,1) AI FF(FIFpop,3)
       FF(FIFpop,1) FF(FIFpop,2) AI
       ];
    % Replace orignal face with 1 (so I can draw them)
    FF(FIFpop,:) = 1;
    J = [J;repmat(FIFpop,3,1)];
    VV = [V;Vpop];

  clf;
  hold on;
  tsurf(FF,VV,'CData',J);
  %tsurf(F,V,'VertexIndices',1,'FaceIndices',1,'FaceColor',orange);
  sct(V(F(FI,1),:).*B(:,1) + V(F(FI,2),:).*B(:,2) + V(F(FI,3),:).*B(:,3),'LineWidth',3);
  sct(VV,'.','SizeData',300);
  %sct(V(uE(uEI,1),:).*uEB(:,1) + V(uE(uEI,2),:).*uEB(:,2),'LineWidth',3,'SizeData',500);
  hold off;
  axis equal;

    % just return? (calls graph coloring many times...)
    return 


    % Remap faces
    break;
  end



  %% Re-express each insertion on unique edge
  %BE = B(IE,:);
  %FIE = FI(IE,:);
  %[~,FJE] = max(BE==0,[],2);
  %uEB = [ ...
  %  BE(sub2ind(size(BE),(1:numel(FJE))',mod(FJE-1+1,3)+1)) ...
  %  BE(sub2ind(size(BE),(1:numel(FJE))',mod(FJE-1+2,3)+1))];
  %E = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  %[uE,IA,EMAP] = unique(sort(E,2),'rows');
  %uEI = EMAP(sub2ind(size(F),FIE,FJE));
  %flip = uE(uEI,1) == E(sub2ind(size(F),FIE,FJE),2);
  %uEB(flip,:) = uEB(flip,[2 1]);
  %% inserting on an edge effects all triangles connected.
  %uE2F = sparse(EMAP,repmat(1:size(F,1),1,3)',1,numel(uE),size(F,1));
  %% Connect all edges in question if they share a face
  %A = (uE2F(uEI,:)*uE2F(uEI,:)')>0;
  %warning('rng');
  %rng(6);
  %C = graph_coloring(A,3);
  %% consider each color
  %for c = reshape(unique(C),1,[])
  %  pop = C==c

  %  % split edges
  %  % remap edge barycentric coordinates
  %  % remap face barycentric coordinates
  %  break;
  %end






end
