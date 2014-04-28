function [S] = readSEL(filename)
  % READSEL read selection from .sel file (as used by Mario Botsch and Olga
  % Sorkine in the accompanying data to their survey "On linear variational
  % surface deformation methods" 
  % http://igl.ethz.ch/projects/deformation-survey/
  %
  % Input:
  %   filename  path to .sel file
  % Output:
  %   S #V list of ids:
  %     0  fixed boundary (silver, completely fixed handle-region)
  %     1  interior (blue, "solved-for", "free" region)
  %     2  handle boundary (gold, re-arrangable handle-region)
  % 
  % Example:
  %  S = readSEL(filename);
  %  % Swap 0 (fixed boundary) with 1 (interior) so that selection mask
  %  % functions as handle ids (with 0 meaning interior)
  %  R = (S<1).*1 + (S==1).*0 + (S>1).*S;
  % 

  fp = fopen(filename,'r');
  % read initial comments
  C = textscan(fp, '#%[^\n]');
  % remaining lines contain one integer per line
  S = fscanf(fp,'%d');
  % clean up file
  fclose(fp);
end
