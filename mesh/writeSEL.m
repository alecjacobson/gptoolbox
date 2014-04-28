function writeSEL(filename,S)
  % writeSEL selection from .sel file (as used by Mario Botsch and Olga
  % Sorkine in the accompanying data to their survey "On linear variational
  % surface deformation methods" 
  % http://igl.ethz.ch/projects/deformation-survey/
  %
  % Input:
  %   filename  path to .sel file
  %   S #V list of ids:
  %     0  fixed boundary (silver, completely fixed handle-region)
  %     1  interior (blue, "solved-for", "free" region)
  %     2  handle boundary (gold, re-arrangable handle-region)
  % 
  % 

  fp = fopen(filename,'w');
  % initial comments
  fprintf(fp,'# per vertex status: 0=fixed, 1=deformable-region, 2=handle\n');
  fprintf(fp,'# %d vertices\n',numel(S));
  % remaining lines contain one integer per line
  fprintf(fp,'%d\n',S(:));
  % clean up file
  fclose(fp);
end

