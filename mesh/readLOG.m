function [V,F,C] = readLOG(filename)
  % READLOG  Read log files output by qslim
  %
  % [V,F,C] = readLOG(filename)
  %
  % Input:
  %   filename  path to .log file
  % Outputs:
  %   V  #V by 3 input mesh vertex positions
  %   F  #F by 3 input mesh triangle indices (1-indexed)
  %   C  #collapse list of edge collapse structs with fields
  %     .a index into V of vertex collapsed to
  %     .b index into V of vertex collapsed from
  %     .da  1 by 3 displacement of vertex a
  %     .rm list of indices in F of faces removed by this collapse
  %     .mod list of indices in F of faces modified by this collapse
  %   
  % See also: perform_edge_collapse

  % first part of file is just the input mesh
  fp = fopen(filename,'r');
  begin = fscanf(fp,' %[^\n]s');
  if(strcmp(begin,'begin')==0)
    error('First line should be "begin"...');
  end
  % force read line feed
  fscanf(fp,'\n');
  V = reshape(fscanf(fp,'v %g %g %g\n'),3,[])';
  F = reshape(fscanf(fp,'f %d %d %d\n'),3,[])';
  endline = fscanf(fp,' %[^\n]s');
  if(strcmp(endline,'end')==0)
    error('Line should be "end"...');
  end
  % force read line feed
  fscanf(fp,'\n');

  c.a =0;
  c.b =0;
  c.rm =[];
  c.mod =[];
  % predict the size of C: reduces time significantly
  C(size(edges(F),1)) = c;
  ci = 1;
  while true
    line = fgets(fp);
    % finished?
    if line==-1
      break;
    end
    [A,count,e,ni] = sscanf(line,'v%% %d %d %g %g %g ',5);
    assert(count==5);
    C(ci).a = A(1);
    C(ci).b = A(2);
    C(ci).da = A(3:5)';
    line = line(ni:end);
    [A,count,e,ni] = sscanf(line,'%d%[^&0-9]');
    % Seems it can be more than 2 removed
    %assert(count==4);
    assert(all(A(2:2:end)==32));
    C(ci).rm = A(1:2:end);
    line = line(ni:end);
    [A,count,e,ni] = sscanf(line,'%c',1);
    assert(count==1);
    assert(A=='&');
    line = line(ni:end);
    [A,count,e,ni] = sscanf(line,'%d');
    C(ci).mod = A(1:end);
    ci = ci+1;
  end
  % only keep real collapses
  C = C(1:ci-1);

  % Now comes the orderedlog of edge collapses
  %
  % This format was backwards engineered.
  %
  % Each line represents an collapse of an edge from vertex b to a
  % v% aid bid dax day daz fa&b1 fa&b2 ... fa&bm & fa|b1 fa|b2 ... fa|bn
  %
  % where aid and bid are the indices of vertices a and b, da* are the
  % coordinates of the displacement of vertex a, fa&b* are the face indices 
  % incident on the edge ab (those removed), and fa|b* are the faces incident
  % on only a or only b (those modified)
  %
  % If F is non-manifold then m may be greater than 2
  %
  % It seems all indices are in terms of the original mesh, as opposed to a
  % cummulative indexing.
  %

  fclose(fp);

end
