function CM = yaron_colormap(n)
  CM = hex2rgb([  
  'ce517e'
  'd16486'
  'd47890'
  'd88a95'
  'dd9c9e'
  'e0ada7'
  'e5beaf'
  'ead1ba'
  'efe2c3'
  'e6e6c4'
  'cfd7bc'
  'b9ccb5'
  'a4c2af'
  '8db7a8'
  '77aca1'
  '54a49e'
  '309c96'
  '159592'
  ]);
  if nargin>=1
    CM = interp1(linspace(0,1,size(CM,1))',CM,linspace(0,1,n)','spline');
  end
end
