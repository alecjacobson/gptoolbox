function [NV,NF,OV,OE,OV1,OE1,epsilon] = remesh_at_handles(V,F,C,P,BE,CE)
  % REMESH_AT_HANDLES  remesh a given mesh using its outline and constraining
  % interrior controls at given handles. Point handles are met exactly and
  % bones are sample at the given rate.
  %
  % [NV,NF] = remesh_at_handles(V,F,C,P,BE,CE)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices
  %  C  list of control vertex positions
  %  P  list of indices into C for point controls, { 1:size(C,1) }
  %  BE  list of bones, pairs of indices into C, connecting control vertices, 
  %    { [] }
  %  CE  list of "cage edges", pairs of indices into ***P***, connecting
  %    control ***points***. A "cage edge" just tells point boundary conditions 
  %    to vary linearly along straight lines between the end points, and to be
  %    zero for all other handles. { [] }
  % Outputs
  %  NV  new list of vertex positions
  %  NF  new list of face indices
  %  OV  list of vertex positions sent to triangle
  %  OE  list of edges sent to triangle
  %  OV1 list of vertex positions before performing close point collapses and
  %    edge splits
  %  OE1 list of edges before performing close point collapses and edge splits
  %  epsilon  used to collapse points and split edges
  %
  % Example:
  %  [NV,NF,OV,OE,OV1,OE1] = remesh_at_handles(V,F,C,P,BE,CE);
  %  % show a plot of what's sent to triangle
  %  subplot(3,1,1);
  %  plot([OV1(OE1(:,1),1) OV1(OE1(:,2),1)]',[OV1(OE1(:,1),2) OV1(OE1(:,2),2)]','-','LineWidth',1);
  %  hold on;
  %  plot(OV1(:,1),OV1(:,2),'.');
  %  hold off;
  %  axis equal;
  %  title('Outline + handle samples');
  %  subplot(3,1,2);
  %  plot([OV(OE(:,1),1) OV(OE(:,2),1)]',[OV(OE(:,1),2) OV(OE(:,2),2)]','-','LineWidth',1);
  %  hold on;
  %  plot(OV(:,1),OV(:,2),'.');
  %  hold off;
  %  axis equal;
  %  title('Input to triangle (collapsed points + split edges)');
  %  % show a plot of the result
  %  subplot(3,1,3);
  %  tsurf(NF,NV);
  %  hold on;
  %  plot([OV(OE(:,1),1) OV(OE(:,2),1)]',[OV(OE(:,1),2) OV(OE(:,2),2)]','-','LineWidth',2);
  %  hold off;
  %  axis equal;
  %  title('Output of triangle');

  %% THIS DOES NOT WORK FOR MESHES THAT AREN"T TOPOLOGICAL DISKS!! HOLES WILL
  %% BE CLOSED

  dim = size(V,2);
  assert(dim == size(C,2));
  % only work with 2D
  assert(dim == 2);

  % Find all edges in mesh, note internal edges are repeated
  E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
  % determine uniqueness of edges
  [u,m,n] = unique(E,'rows');
  % determine counts for each unique edge
  counts = accumarray(n(:), 1);
  % extract edges that only occurred once
  O = u(counts==1,:);
  % extract unique vertex indices on outline
  [u,m,n] = unique(O(:));
  % original map O = IM(O)
  IM = 1:size(V,1);
  IM(O(:)) = n;
  OV = V(u,:);
  OE = IM(O);



  samples_per_edge = 10;
  % Bone samples
  [BES,BESE] = sample_edges(C,BE,samples_per_edge);
  % Cage edge samples
  [CES,CESE] = sample_edges(C(P,:),CE,samples_per_edge);
  % Append handle samples
  OV = [OV;C(P,:);BES;CES];
  OE = [ ...
    OE; ...
    (size(OV,1)-size(BES,1)-size(CES,1))+BESE; ...
    (size(OV,1)-size(CES,1))+CESE];

  % save old values
  OV1 = OV;
  OE1 = OE;

  % phony counts to enter while loop
  prev_v_count = size(OV,1)+1;
  prev_e_count = size(OE,1)+1;
  % this epsilon is a parameter that probably needs to be tweaked and should at
  % least be exposed.
  % Use fraction of average edge length
  epsilon = mean(sqrt(sum((OV(OE(:,1),:) - OV(OE(:,2),:)).^2,2)))/8;
  % continue to collapse points and split edges until nothing more to do
  while( prev_v_count ~= size(OV,1) || prev_e_count ~= size(OE,1))
    prev_v_count = size(OV,1);
    prev_e_count = size(OE,1);
    [OV,OE] = collapse_close_points(OV,OE,epsilon);
    [OV,OE] = snap_points_to_close_edges(OV,OE,epsilon);
  end

  % Again, max area term here is a heuristic
  %sum((OV(OE(:,1)) - OV(OE(:,2))).^2,2)
  % use a multiple of min edge length
  %max_area = 2*(sqrt(3)/4)*min(sum((OV(OE(:,1),:) - OV(OE(:,2),:)).^2,2));
  max_area = 8*(sqrt(3)/4)*min(sum((OV(OE(:,1),:) - OV(OE(:,2),:)).^2,2));
  [NV,NF] = triangle(OV,OE,[],'Quality',30,'MaxArea',max_area);
end
