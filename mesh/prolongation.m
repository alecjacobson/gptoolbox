function [P] = prolongation(CV,CF,V,varargin)
  % PROLONGATION Build a linear prolongation operator taking solutions on a
  % coarse mesh (CV,CF) and prolongating (upsampling) them onto vertices of a
  % fine mesh (V). Useful for multigrid and embedded mesh deformations.
  %
  % P = prolongation(CV,CF,V)
  %
  % Inputs:
  %   CV  #CV by dim list of coarse mesh vertex positions
  %   CF  #CF by dim+1 list of coarse mesh element indices into CV
  %   V  #V by dim list of fine mesh vertex positions
  %   Optional:
  %     'Extrapolation'  followed by name of method to use for points in V
  %       lying outside (CV,CF). One of the following:
  %         {'linear'}   finds closest element and uses barycentric coordinates
  %           (with some negative weights). This is linearly precise.
  %         'contstant'  finds closest point and uses its barycentric
  %           coordinates (all non-negtive since it lies on a facet). This is
  %           only constantly precise.
  % Outputs:
  %   P  #V by #CV prolongation matrix so that X = P * CX prolongates a
  %     solution CX on the coarse mesh to a solution X on the fine mesh.
  % 

  [I] = in_element_aabb(CV,CF,V);
  I0 = find(I==0);
  extrapolation = 'linear';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Extrapolation'}, ...
    {'extrapolation'});
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

  switch size(CF,2)
  case 4
    [BF,BJ,BK] = boundary_faces(CF);
    BV = V;
    % Extrapolate
    % Snap to boudnary and use snapped points interpolation
    [~,BI,BVI0] = point_mesh_squared_distance(V(I0,:),CV,BF);
    I(I0) = BJ(BI);
    switch extrapolation
    case 'linear'
      % Leave BV
    case 'constant'
      BV(I0,:) = BVI0;
    otherwise
      error('Unknown extrapolation (%s)',extrapolation);
    end
    [IV] = max(I,[],2);
    % recover barycentric coordinates (interpolation & extrapolation)
    B = barycentric_coordinates( ...
      BV,CV(CF(IV,1),:),CV(CF(IV,2),:),CV(CF(IV,3),:),CV(CF(IV,4),:));
    %% Snap to closest vertex
    %B(any(B<0,2),:) = bsxfun(@eq,max(B(any(B<0,2),:),[],2),B(any(B<0,2),:));
  case 3
    % Dirty trick to make V,CV 3D
    BV = V;
    if ~isempty(I0)
      pV = V;
      pCV = CV;
      pBV = BV;
      pV(:,end+1:3) = 0;
      pCV(:,end+1:3) = 0;
      [~,I(I0),pBVI0] = point_mesh_squared_distance(pV(I0,:),pCV,CF);
      switch extrapolation
      case 'linear'
        % Leave BV
      case 'constant'
        BV(I0,:) = pBVI0(:,1:2);
      otherwise
        error('Unknown extrapolation (%s)',extrapolation);
      end
    end
    [IV] = max(I,[],2);
    B = barycentric_coordinates( ...
      BV,CV(CF(IV,1),:),CV(CF(IV,2),:),CV(CF(IV,3),:));
  end

  P = sparse( ...
    repmat(1:size(V,1),size(CF,2),1)', ...
    CF(IV,:), ...
    B, ...
    size(V,1),size(CV,1));
end
