function writeMDD(mdd_filename,VV)
  % WRITEMDD Write a mesh animation to a .mdd file, e.g., for use in Blender
  %
  % Inputs:
  %   mdd_filename  path to .mdd file
  %   VV   #V by 3 by #frames mesh animation
  % 
  fp = fopen(mdd_filename,'w','ieee-be');
  fwrite(fp,size(VV,3),'uint32');
  fwrite(fp,size(VV,1),'uint32');
  fwrite(fp,linspace(0,1,size(VV,3)),'float');
  fwrite(fp,permute(VV,[2 1 3]),'float');
  fclose(fp);
end
