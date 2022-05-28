function [mV,mF,VV,FF] = meshfix_components(V,F)
  % MESHFIX_COMPONENTS
  %
  % [mV,mF,VV,FF] = meshfix_components(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of surface mesh vertex positions
  %   F  #F by 3 list of surface mesh face indices into V
  % Outputs:
  %   mV  #mV by 3 list of output surface mesh vertex positsions
  %   mF  #mF by 3 list of output surface mesh face indices into SV
  %

  C = connected_components(F); C = C(F(:,1));
  k = max(C);
  VV = cell(k,1);
  FF = cell(k,1);
  mV = cell2mat(VV);
  for i = 1:k
    Fi = F(C==i,:);
    [Vi,I] = remove_unreferenced(V,Fi);
    Fi = reshape(I(Fi),size(Fi));
    try
      [VV{i},FF{i}] = meshfix(Vi,Fi);
    catch
      VV{i} = zeros(0,3);
      FF{i} = zeros(0,3);
    end
  end
  mV = cell2mat(VV);
  n = cumsum([0 cellfun(@(V) size(V,1),{VV{:}})]);
  mF = cell2mat(arrayfun(@(i) n(i)+FF{i},1:numel(VV),'UniformOutput',false)');
end
