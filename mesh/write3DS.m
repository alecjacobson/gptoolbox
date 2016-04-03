function write3DS(filename,V,F)
  % WRITE3DS Write a mesh (V,F) to a .3ds file
  %
  % write3DS(filename,V,F)
  % 
  % Inputs:
  %   filename  path to .3ds file
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into V
  % 
  % Warning: outputs files are no readable by meshlab
  %
  n = size(V,1);
  m = size(F,1);

  % sizes of chunks below
  f_len = 10+2*m*4;
  v_len = 8+4*n*3;
  t_len = v_len+f_len+6;
  name = ['mesh' 0];
  name_len = numel(name);
  m_len = name_len+t_len+6;
  s_len = 10;
  mver_len = 10;
  e_len = s_len + mver_len + m_len+6;
  ver_len = 10;
  h_len = 2+6;
  on_len = 11+6;
  j_len = 18+20+20+20;
  mi_len = h_len + j_len + on_len + 6;
  k_len = mi_len + 6;
  main_len = k_len + ver_len + e_len+6;

  f = fopen(filename, 'wb');
  vertices_id = 16656;
  texture_coords_id = 16704;
  faces_id = 16672;
  % main
  fwrite(f,19789,'ushort');
  fwrite(f,main_len,'int');
  % version
  fwrite(f,2,'ushort');
  fwrite(f,ver_len,'int');
  fwrite(f,3,'int');
  % 3d editor
  fwrite(f,15677,'ushort');
  fwrite(f,e_len,'int');
  % version
  fwrite(f,15678,'ushort');
  fwrite(f,mver_len,'int');
  fwrite(f,3,'int');
  % scale?
  fwrite(f,256,'ushort');
  fwrite(f,s_len,'int');
  fwrite(f,1.0,'float32');
  % mesh
  fwrite(f,16384,'ushort');
  fwrite(f,m_len,'int');
  fwrite(f,name,'char');
  % triangle mesh
  fwrite(f,16640,'ushort');
  fwrite(f,t_len,'int');
  % vertices
  fwrite(f,16656,'ushort');
  fwrite(f,v_len,'int');
  fwrite(f,n,'ushort');
  fwrite(f,V','float32');

  % faces
  F = [F-1 zeros(size(F,1),1)];
  fwrite(f,16672,'ushort');
  fwrite(f,f_len,'int');
  fwrite(f,m,'ushort');
  fwrite(f,F','ushort');
  fwrite(f,0,'ushort');

  % For some reason meshlab refuses to read the .3ds file unless it has this
  % junk at the end

  % keyframer
  fwrite(f,45056,'ushort');
  fwrite(f,k_len,'int');
  % mesh information
  fwrite(f,45058,'ushort');
  fwrite(f,mi_len,'int');
  % hiearchy position
  fwrite(f,45104,'ushort');
  fwrite(f,h_len,'int');
  fwrite(f,0,'ushort');
  % object name
  fwrite(f,45072,'ushort');
  fwrite(f,on_len,'int');
  fwrite(f,['mesh' 0],'char');
  fwrite(f,0,'int');
  fwrite(f,-1,'short');
  % Junk
  % Object Pivot Point
  fwrite(f,45075,'ushort');
  fwrite(f,12+6,'int');
  fwrite(f,zeros(12,1),'char');
  % position track
  fwrite(f,45088,'ushort');
  fwrite(f,14+6,'int');
  fwrite(f,zeros(14,1),'char');
  % rotation track
  fwrite(f,45089,'ushort');
  fwrite(f,14+6,'int');
  fwrite(f,zeros(14,1),'char');
  fwrite(f,45090,'ushort');
  fwrite(f,14+6,'int');
  fwrite(f,zeros(14,1),'char');
  

end
