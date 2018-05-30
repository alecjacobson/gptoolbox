function [V,T,F] = readMSH(filename)
  % Input:
  %   filename path to gmsh .msh file
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   T  #T by 3 list of tet mesh indices into V
   
  %function expect_line(f,str)
  %  line = fscanf(f,'%s\n',1);
  %  if strcmp(line,str) ~= 1
  %    assert(false,'"%s" ~= "%s"',line,str);
  %  end
  %end

  %f = fopen(filename,'r');
  %expect_line(f,'$MeshFormat');
  %vfd = fscanf(f,'%g %g %g\n',3);
  %version = vfd(1);
  %filetype = vfd(2);
  %datasize = vfd(3);
  %assert(datasize == 8,'only doubles are supported');
  %ascii = filetype == 0;
  %if ascii
  %  assert(false,'Not implemented');
  %else
  %  one = fread(f,1,'*int');
  %  assert(one == 1,'Unexpected endian, not supported');
  %  expect_line(f,'$EndMeshFormat');
  %  expect_line(f,'$Nodes');
  %  n = fscanf(f,'%d\n',1);
  %  id1 = fread(f,1,'*int');
  %  V = fread(f,[3 n],'3*double',4);
  %  fseek(f,-4,0);
  %  expect_line(f,'$EndNodes');
  %  expect_line(f,'$Elements');
  %  num_elements_total = fscanf(f,'%d\n',1);
  %  header = fread(f,3,'*int');
  %  ss = header(1)
  %  m = header(2)
  %  num_tags = header(3)
  %  id1 = fread(f,1,'*int');
  %end
  %fclose(f);
end
