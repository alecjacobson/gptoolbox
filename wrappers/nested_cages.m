function [ncV,ncF] = nested_cages(V,F,dV,dF)
  % [ncV,ncF] = nested_cages(V,F,dV,dF)
  %
  % 
  % Example:
  %   [V,F] = subdivided_sphere(2);
  %   [dV,dF] = subdivided_sphere(0);
  %   dV=dV*1.2;clf;
  %   [ncV,ncF] = nested_cages(V,F,dV,dF);
  %
  path_to_nested_cages = '/Users/alecjacobson/Repos/nested_cages/build/nested_cages';

  tmp = tempname;
  input_filename = [tmp 'nested_cages_input.ply'];
  start_filename = [tmp 'nested_cages_start.ply'];
  writePLY(input_filename,V,F);
  writePLY(start_filename,dV,dF);
  prefix = [tmp 'nested_cages_output'];
  cmd = sprintf('%s %s 2 %s None Volume %s',path_to_nested_cages,input_filename,start_filename,prefix);
  [status,result] = system(cmd);
  if status ~= 0
    warning(result)
  end
  [ncV,ncF] = load_mesh([prefix '_1.obj']);
  delete([prefix '_1.obj'])
  delete([prefix '_0.obj'])
  delete(input_filename);
  delete(start_filename);
end

