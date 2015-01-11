function [VV,FF] = extrude(V,F,varargin)
  % EXTRUDE Extrude a 2d mesh in the z direction by 1, connecting boundaries
  % apropriately 
  %
  % [VV,FF] = extrude(V,F)
  % [VV,FF] = extrude(V,F,'ParmeterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 2 list of 2d vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'Cap' followed by {true} of false whether to put a *bottom cap*
  % Outputs:
  %   VV #VV by 3 list of 3d vertex positions
  %   FF  #FF by 3 list of triangle indices into VV
  %

  cap = true;
  levels = 1;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Cap','Levels'}, ...
    {'cap','levels'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace 
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  assert(size(V,2) == 2, 'Vertices should be 2d');
  %% Copy top and bottom
  %  VV = [V 1+0*V(:,1);V 0*V(:,1)];
  switch size(F,2)
  case {3,4}
    % connect boundaries
    O = outline(F);
  case 2
    assert(~cap,'Cap should be false when input is curve');
    % Curve as input
    FF = [];
    O = F;
  end
  assert(size(unique(sort(O,2),'rows'),1) == size(O,1));

  % Rearrange vertices on bottom
  [VB,FB,OB,IM] = faces_first(V,F,O);
  n = size(V,1);

  no = numel(unique(OB));

  VL = [ ...
      repmat(VB(1:max(OB(:)),:),[levels 1]) ...
        reshape(repmat(linspace(1-1/levels,0,levels),[no,1]),[no*levels 1])];
  if cap
    VV = [VB ones(n,1);VL(1:no*(levels-1),:);VB zeros(n,1)];
  else
    VV = [VB ones(n,1);VL];
  end
  FO = ones(2*size(OB,1)*levels,size(F,2));
  for level = 1:levels
    if level == 1
      off_t = 0;
    else
      off_t = n+(level-2)*no;
    end
    off_b = n+(level-1)*no;
    switch size(F,2)
    case 4
      FO((level-1)*size(OB,1)+(1:size(OB,1)),:) = [ ...
        off_t+OB(:,[2 1]) off_b+OB(:,1:2)];
    case 3
      FO((level-1)*2*size(OB,1)+(1:size(OB,1)*2),:) = [ ...
        off_b+OB(:,1) off_t+OB(:,[2 1]); ...
        off_t+OB(:,2) off_b+OB(:,1:2)];
    end
  end
  
  % remap so that faces on planes are consistent with input
  if cap
    FF = [FB;FO;n+(levels-1)*no+fliplr(FB)];
    RIM = 1:(n+(levels-1)*no+n);
    RIM(IM) = 1:n;
    RIM(n+IM) = n+(1:n);
  else
    FF = [FB;FO];
    RIM = 1:(n+levels*no);
    RIM(IM) = 1:n;
  end
  VV(RIM,:) = VV;
  FF = RIM(FF);


end
