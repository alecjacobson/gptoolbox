function [iter,dt,t,x,R,v,omega,M,I,fixed] = readSCISIM(h5_filename)

  path_to_print_meshes = ...
    '/Users/ajx/Dropbox/scisim/post_processing/print_meshes.py';
  path_to_python = '/usr/local/bin/python';
  out_filename = tempname;
  cmd = sprintf('%s %s %s >%s', ...
    path_to_python,path_to_print_meshes,h5_filename,out_filename);
  [status,r] = system(cmd);
  if status ~= 0
    error(sprintf('\n%s\n\n%s',cmd,r));
  end

  fp = fopen(out_filename,'r');
  hash = fscanf(fp,'Git hash: %s\n',1);
  iter = fscanf(fp,'Iteration: %d\n',1);
  dt = fscanf(fp,'Timestep: %f\n',1);
  t = fscanf(fp,'Time: %f\n',1);

  x     = zeros(100,3);
  R     = zeros(100,9);
  v     = zeros(100,3);
  omega = zeros(100,3);
  M     = zeros(100,1);
  I     = zeros(100,9);
  fixed = zeros(100,1);
  mesh_name = cell(1,100);

  count = 0;
  while true
    id = fscanf(fp,'Body: %d\n',1);
    if isempty(id)
      break;
    else 
      count = count+1;
      assert(count == id+1);
    end
    x(count,:) = fscanf(fp,'  x: [ %f %f %f ]\n');
    R(count,1:3) = fscanf(fp,'  R: [%f %f %f]\n',3);
    R(count,4:6) = fscanf(fp,'     [%f %f %f]\n',3);
    R(count,7:9) = fscanf(fp,'     [%f %f %f]\n',3);
    v(count,:) = fscanf(fp,'  v: [ %f %f %f ]\n');
    omega(count,:) = fscanf(fp,'  omega: [ %f %f %f ]\n');
    M(count) = fscanf(fp,'  M: %f\n',1);
    I(count,1:3) = fscanf(fp,'  I: [%f %f %f]\n',3);
    I(count,4:6) = fscanf(fp,'     [%f %f %f]\n',3);
    I(count,7:9) = fscanf(fp,'     [%f %f %f]\n',3);
    fixed(count) = fscanf(fp,'  fixed: %f\n',1);
    mesh_name{count} = fscanf(fp,'  mesh_name: %s\n',1);

    if count==size(x,1)
      x     = [x    ;zeros(count,3)];
      R     = [R    ;zeros(count,9)];
      v     = [v    ;zeros(count,3)];
      omega = [omega;zeros(count,3)];
      M     = [M    ;zeros(count,1)];
      I     = [I    ;zeros(count,9)];
      fixed = [fixed;zeros(count,1)];
      mesh_name = {mesh_name{:} cell(1,count)};
    end
  end

  x = x(1:count,:);
  R = R(1:count,:);
  v = v(1:count,:);
  omega = omega(1:count,:);
  M = M(1:count,:);
  I = I(1:count,:);
  fixed = fixed(1:count,:);
  fclose(fp);
  delete(out_filename);
end
