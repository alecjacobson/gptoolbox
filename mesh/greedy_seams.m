function Ebest = greedy_seams(V,F)
  % Ebest = greedy_seams(V,F)
  %
  % An attempt to cut a mesh into seams
  %
  % Inputs:
  %  V  #V by 3 list of vertex positions
  %  F  #F by 3|4 list of triangle or quad indices
  % Outputs:
  %  Ebest   #E by 2 list of edges in the best seam set found

  % handle quads, too
  if size(F,2) == 4
    F = [Q(:,[1 2 3]);Q(:,[1 3 4])];
    allE = [Q(:,[1 2]);Q(:,[2 3]);Q(:,[3 4]);Q(:,[4 1])];
  else
    allE  = [F(:,2:3);F(:,[3 1]);F(:,1:2)];
  end

  V = (V-min(V))/max(max(V)-min(V));
  rng(0);
  [~,I] = farthest_points(V,5,'Distance','biharmonic','F',F);
  
  E = unique(sort(reshape(allE,[],2),2),'rows');
  [A,AE] = adjacency_incident_angle_matrix(V,E);
  [AI,AJ] = find(AE>0);
  AH = A(sub2ind(size(A),AI,AJ));
  l = edge_lengths(V,E);
  % coming into or leaving an edge pay ½ cost of length
  AL = (l(AI)+l(AJ))*0.5;
  % rebuild A
  AC = AL+2*mean(AL)*(-log(AH)+log(pi));
  A = sparse(AI,AJ,AC,size(E,1),size(E,1));
  % special directed graph
  V2E = sparse(E,repmat(1:size(E,1),2,1)',[l l],size(V,1),size(E,1));
  n = size(V,1);
  k = size(E,1);
  % source,sink,edges
  A = [sparse(n,2*n) V2E;sparse(n,2*n+k);sparse(k,n) V2E' A];
  % special super-sink
  A(end+1,end+1) = 0;
  
  %    clf;
  %    hold on;
  %    tsh = tsurf(F,V,'FaceVertexCData',repmat(blue,size(V,1),1),'EdgeColor',blue*0.3,falpha(0.7,1.0),fphong,fsoft);
  %    orange = [0.9 0.6 0.0];
  %    esh = tsurf([1 1 1],V+1e-3*per_vertex_normals(V,F),'LineWidth',3,'EdgeColor',orange);
  %    hold off;
  %    axis equal;
  %    camproj('persp');
  %    set(gca,'Position',[0 0 1 1],'Visible','off');set(gcf,'Color','w');
  %    view(-51,50);
  %    l = add_lights();
  %    add_shadow(tsh,l{5});
  %    %apply_ambient_occlusion(tsh,'SoftLighting',false,'AddLights',false);
  %    hold on;
  %    ssh = sct(V(I,:),'o','SizeData',100,'MarkerFaceColor',orange,'LineWidth',2,'MarkerEdgeColor',orange*0.5);
  %    hold off;
  %    view(-51,18);
  %    drawnow;
  
  P = perms(1:size(I,1));
  rng(0);
  P = P(randperm(end),:);
  size(P)
  min_cost = inf;
  Ebest = [];
  for p = 1:size(P,1)
    progressbar(p,size(P,1));
    Ip = I(P(p,:));
    %Ip = I([4 3 1 5 2]);
  
    Erun = zeros(0,2);
    S = [];
    % clear super-sink
    A(:,end) = 0;
    cost = 0;
    for i = reshape(Ip,1,[])
      if isempty(S)
        S = i;
      else
        % shortest path to super-sink
        [~,path] = graphshortestpath(A,i,size(A,1));
        cost = cost + sum(A(sub2ind(size(A),path(1:end-1),path(2:end))));
        Ei = E(path(2:end-2)-2*n,:);
        %if cost > min_cost
        %  % abandon early.
        %  break;
        %end
        Erun = [Erun;Ei];
        S = union(S,Ei(:));
      end
      A(n+S,end) = 1;
    end
    %esh.Faces = Erun(:,[1 2 2]);
    %drawnow;
    if cost < min_cost
      min_cost = cost
      Ebest = Erun;
    end
  end
  
  esh.Faces = Ebest(:,[1 2 2]);
  drawnow;
  
  [F,I] = cut_edges(F,Ebest);
  V = V(I,:);
end
