function [VV,FF] = extrude(V,F,varargin)
  % EXTRUDE Extrude a 2d mesh in the z direction by 1, connecting boundaries apropriately 
  %
  % [VV,FF] = extrude(V,F)
  % [VV,FF] = extrude(V,F,'ParmeterName',ParameterValue,...)
  %
  % Inputs:
  %  V  #V by 2 list of 2d vertex positions
  %  F  #F by 3 list of triangle indices into V
  %  Optional:
  %    'Cap' followed by {true} of false whether to put a *bottom cap*
  % Outputs:
  %  VV #VV by 3 list of 3d vertex positions
  %  FF  #FF by 3 list of triangle indices into VV
  %

  cap = true;
  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'Cap'
      v = v+1;
      assert(v<=numel(varargin));
      cap = varargin{v};
    end
    v = v+1;
  end

  % Copy top and bottom
  VV = [V 1+0*V(:,1);V 0*V(:,1)];
  switch size(F,2)
  case 3
    % Mesh as input
    if cap
      FF = [F;size(V,1)+fliplr(F)];
    else
      FF = F;
    end

    % connect boundaries
    O = outline(F);
  case 2
    % Curve as input
    FF = [];
    O = F;
  end
  assert(size(unique(sort(O,2),'rows'),1) == size(O,1));

  FO = [ ... 
    size(V,1)+O(:,1) O(:,[2 1]); ...
    O(:,2) size(V,1)+O(:,1:2)];

  FF = [FF;FO];

end
