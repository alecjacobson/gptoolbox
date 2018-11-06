function [VV,FF,RV,RF,PRV,PE] = joints(V,E,varargin)
  % JOINTS Construct joints around a wire mesh given as a graph
  % 
  % [VV,FF,RV,RF,PRV,PE] = joints(V,E);
  % [VV,FF,RV,RF,PRV,PE] = joints(V,E,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of wire graph vertex positions
  %   E  #E by 2 list of edge indices into V
  %   Optional:
  %     'Radius'  followed by dowel rod radius
  %     'Tol'  followed by engineering tolerance
  %     'PolySize'  followed by number of sides on dowel rod cross-sectional
  %       polygon
  %     'JointThickness' followed by thickness of joints {0.25*Radius}
  % Outputs:
  %   VV  #VV by 3 list of vertex positions of joint mesh
  %   FF  #FF by 3 list of triangle indices into VV of joint mesh
  %   VV  #VV by 3 list of vertex positions of rod mesh
  %   FF  #FF by 3 list of triangle indices into VV of rod mesh
  %   PRV  2*#E by 3  list of offset edge endpoint positions
  %   PE   #E by 2 list of edge indices into PRV
  %   

  % default values
  poly = 5;
  r = 0.05;
  tol = 0.01;
  th = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'JointThickness','PolySize','Radius','Tol'}, ...
    {'th','poly','r','tol'});
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

  if isempty(th)
    th = 0.25*r;
  end

  % Create wire mesh used for joints
  fprintf('wiremesh...\n');
  [WV,WF,WJ] = wire_mesh(V,E,'PolySize',poly,'Thickness',2*(r+th+tol));
  
  % Edge faces
  EF = WF(WJ>size(V,1),:);
  EJ = WJ(WJ>size(V,1),:)-size(V,1);
  % Edge Vertices
  EV = WV(EF(:),:);
  EI = repmat(EJ,3,1);
  
  % Edge root
  ER = V(E(:,1),:);
  % unit edge vector
  EU = normalizerow(V(E(:,2),:)-V(E(:,1),:));
  % project each onto its respect edge vector
  ED = sparse(EI,(1:size(EV,1))',sum((EV-ER(EI,:)).*EU(EI,:),2),size(E,1),size(EV,1));
  % Amount of overhang
  hang = 2*r;
  Emin = minnz(ED')';
  Emax = max(ED,[],2);
  
  PRV = [ER + Emin.*EU;ER+Emax.*EU];
  PIV = [ER + (Emin-th).*EU;ER+(Emax+th).*EU];
  PV = [ER + (Emin+hang).*EU;ER+(Emax-hang).*EU];
  PE = reshape(1:2*size(E,1),size(E,1),2);

  % Over hang
  % https://math.stackexchange.com/a/1001549/35376
  [CV,CF,CJ] = edge_cylinders( ...
    PV,PE, ...
    'PolySize',poly, ...
    'Thickness',2*(1+2*(1-cos(pi/poly)))*(r+th+tol));

  % inner bar
  [IV,IF] = edge_cylinders(PIV,PE, 'PolySize',poly, 'Thickness',2*(r+tol));
  fprintf('booleaning...\n');
  [VV,FF] = mesh_boolean(WV,WF,[CV;IV],[CF;size(CV,1)+IF],'minus');
  [RV,RF] = edge_cylinders(PRV,PE, 'PolySize',poly, 'Thickness',2*r);

end

