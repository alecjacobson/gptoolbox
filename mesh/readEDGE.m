function [E,RP,RD] = readEDGE(filename)
  % READEDGE read .edge file used by Triangle program to output edges of
  % voronoi diagrams
  %
  % [E,RP,RD] = readedge(filename)
  % 
  % Input:
  %   filename  input file name
  %
  % Output:
  %   E  #E by 2 list of voronoi internal edges
  %   RP  #RP by 1 list of voronoi infinite ray start points
  %   RD  #RP by 2 list of voronoi infinite direction vectors
  %
  % Example:
  %   D = readNODE(data_node_file);
  %   V = readNODE(voronoi_node_file);
  %   [E,RP,RD] = readEDGE(voronoi_edge_file);
  %   RD = normalizerow(RD);
  %   offset = 2*max(sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2)));
  %   EX = [V(E(:,1),1),V(E(:,2),1); V(RP,1) V(RP,1)+offset*RD(:,1)]';
  %   EY = [V(E(:,1),2),V(E(:,2),2); V(RP,2) V(RP,2)+offset*RD(:,2)]';
  %   h = plot(EX,EY,'-',D(:,1),D(:,2),'.');
  %   set(h(1:end-1),'xliminclude','off','yliminclude','off');
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also readNODE
  %

  edge_file_handle = fopen(filename);
  header = fscanf(edge_file_handle,'%d %d',2);
  % number of edge lines
  num_e = header(1);
  % number of attributes
  num_a = header(2);
  if(num_a ~= 0)
    error('Attributes are not supported.');
  end

  E = [];
  RP = [];
  RD = [];
  min_e = inf;
  while(true)
    % read next whole line
    line = fscanf(edge_file_handle,' %[^\n]s');
    [e,count] = sscanf(line,'%d %d %d %g %g %g',6);
    if(count == 3)
      E = [E;e(2) e(3)];
    elseif(count == 6)
      assert(e(3) == -1);
      RP = [RP; e(2)];
      RD = [RD; e(4) e(5) e(6)];
    elseif(count == 5)
      assert(e(3) == -1);
      RP = [RP; e(2)];
      RD = [RD; e(4) e(5)];
    else
      break;
    end
    min_e = min([e(1) min_e]);
  end
  % adjust for 0-based indices
  if(min_e == 0)
    E = E + 1;
    RP = RP + 1;
  end
  fclose(edge_file_handle);
end

