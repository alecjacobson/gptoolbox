function [f,G,H] = symmetric_dirichlet(U,F,P,varargin)
  % SYMMETRIC_DIRICHLET Compute the symmetric Dirichlet energy (and its
  % derivatives) for 2D parameterization of a given *strictly feasible** 2D
  % mesh.  See "Bijective Parameterization with Free Boundaries" [Smith and
  % Schaefer 2015].
  %
  % Inputs:
  %   U  #V by 2 list of input 2D mesh vertex displacements
  %   F  #F by 3 list of triangle indices into rows of U/P
  %   P  #P by 2|3 list of input 3D mesh vertex positions
  %    or
  %      #F*3 by 2|3 list of input 3D mesh corner positions: P = V(F,:)
  %      e.g., the output of plane_project.m
  % Outputs:
  %   f  scalar total objective value
  %   G  #U*2 gradient 
  %   H  #U*2 by #U*2 sparse Hessian matrix
  %
  % Note: this function uses a special symbolic library trick which generates
  % file symmetric_dirichlet_sym.m if it doesn't already exist. If your change
  % the symbolic-math part of this file, you must delete this automatically
  % generated file.
  %

  psd_project = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PSD'}, ...
    {'psd_project'});
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


  % Add column of zeros if necessary
  P(:,end+1:3) = 0;
  if size(P,1) == size(U,1)
    % [ x x x  y y y z z z]
    PI = reshape(F(:)+(0:size(P,2)-1)*size(P,1),size(F,1),size(P,2)*size(F,2));
    %P = reshape(P(F,:),size(F,1),[]);
    P = P(PI);
  elseif size(P,1) == numel(F)
    P = reshape(P,size(F,1),[]);
  end

  % [ x x x  y y y]
  nu = numel(U);
  UI = reshape(F(:)+(0:size(U,2)-1)*size(U,1),size(F,1),size(U,2)*size(F,2));
  U = U(UI);

  % Writing the symbol
  path = mfilename('fullpath');
  aux = [path '_sym.m'];
  % should also check date...
  if ~exist(aux,'file')
    % From here, only touch X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sU = sym('U',[1 6]); % [x₁ x₂ x₃ y₁ y₂ y₃]
    sP = sym('P',[1 9]); % [x₁ x₂ x₃ y₁ y₂ y₃ z₁ z₂ z₃]
    % [sU1,sU2,sU3] = deal(sU(1:3:end),sU(2:3:end),sU(3:3:end));
    % [sP1,sP2,sP3] = deal(sP(1:3:end),sP(2:3:end),sP(3:3:end));
    sU1 = sU(1:3:end); % [u₁ u₄] = [x₁ y₁]
    sU2 = sU(2:3:end); % [u₂ u₅] = [x₂ y₂]
    sU3 = sU(3:3:end); % [u₃ u₆] = [x₃ y₃]
    sP1 = sP(1:3:end); % [p₁ p₄ p₇] = [x₁ y₁ z₁]
    sP2 = sP(2:3:end);
    sP3 = sP(3:3:end);
    
    % squared area
    sD2U = (0.5*((sU2-sU1)*[0 1;-1 0]*(sU3-sU1).')).^2;
    sqrlen = @(sX) sum(sX.^2,2);
    sD2P = sqrlen(0.5*cross(sP2-sP1,sP3-sP1,2));
    sf = (1+sD2P./sD2U) .*  ( ...
      (sqrlen(sU3-sU1).*sqrlen(sP2-sP1) + sqrlen(sU2-sU1).*sqrlen(sP3-sP1))./ ...
        (4*sqrt(sD2P)) - ...
      (sum((sU3-sU1).*(sU2-sU1),2).*sum((sP3-sP1).*(sP2-sP1),2))./ ...
        (2*sqrt(sD2P)));
  
    hess = @(sf,sX) cell2sym(arrayfun(@(g) gradient(g,sX),gradient(sf,sX),'UniformOutput',false));
  
    aux_handle = ...
      matlabFunction(sf,gradient(sf,sU),hess(sf,sU),'vars',{sU,sP},'File',aux);
  else
    aux_name = [mfilename('func') '_sym'];
    aux_handle = str2func(aux_name);
  end

  faux=[];gaux = [];Haux = [];
  switch nargout
  case {0,1}
    [faux] = aux_handle(U,P);
  case 2
    [faux,gaux] = aux_handle(U,P);
  case 3
    [faux,gaux,Haux] = aux_handle(U,P);
  end

  % unnecessary indirection?
  f_fun = @(U) faux;
  dfdU_fun = @(U) gaux;
  d2fdU2_fun = @(U) Haux;
  f = sum(f_fun(U));
  if nargout<=1
    return;
  end
  dfdU = dfdU_fun(U);
  G = full_sparse(UI,ones(size(UI)),reshape(dfdU,size(UI)),nu,1);
  if nargout<=2
    return;
  end
  d2fdU2 = double(d2fdU2_fun(U));

  d2fdU2 = reshape(d2fdU2,[],6*6);
  if psd_project
    d2fdU2 = psd_project_rows(d2fdU2);
  end

  HI = repmat(UI,[1 1 size(UI,2)]);
  HJ = permute(repmat(UI,[1 1 size(UI,2)]),[1 3 2]);
  H = fast_sparse(HI(:),HJ(:),d2fdU2(:),nu,nu);



end
