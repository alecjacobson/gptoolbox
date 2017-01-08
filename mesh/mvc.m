function [W,TA] = mvc(V,P,E,varargin)
  % MVC Compute mean value coordinates given a list of query points, list of
  % corner locations and list of oriented edges indexing those corners.
  %
  % [W,TA] = mvc(V,P,E,varargin)
  %
  % Inputs:
  %   V  #V by 2 list of query points
  %   P  #P by 2 list of corner positions
  %   E  #E by 2 list of oriented edge indices into P
  %   Optional:
  %     'EnforceBoundaryConditions'  followed by whether to use special
  %       handling for query points in V that lie on or near the control
  %       segments. This can be slow, so if you're sure that V are far from
  %       (P,E) then you can safely set this to false. {true}
  % Outputs:
  %   W  #V by #P list of mean value weights
  %  
  %

  enforce_boundary_conditions = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'EnforceBoundaryConditions'}, ...
    {'enforce_boundary_conditions'});
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

  assert(size(V,2) == size(P,2),'V and P need to be in same dimension');
  assert(size(V,2) == 2);
  % Vectors from queries to cage vertices
  %VP = bsxfun(@minus, permute(P,[3 1 2]), permute(V,[1 3 2]));
  %R = sqrt(sum(VP.^2,3));
  % Faster to avoid 3D arrays
  VP1 = bsxfun(@minus,P(:,1)',V(:,1));
  VP2 = bsxfun(@minus,P(:,2)',V(:,2));
  R = sqrt(VP1.^2 + VP2.^2);
  %% sine of alpha for each edge
  %SA = (VP1(:,E(:,1)).*VP2(:,E(:,2)) -VP1(:,E(:,2)).*VP2(:,E(:,1)))./ ...
  %  (R(:,E(:,2)).*R(:,E(:,1)));
  %% cosine of alpha for each edge
  %CA = sum(VP(:,E(:,2),:).* VP(:,E(:,1),:),3) ./ (R(:,E(:,2)).*R(:,E(:,1)));
  %% https://en.wikipedia.org/wiki/Tangent_half-angle_formula
  %TA = SA./(1+CA);
  % More accurate
  %TA = tan((atan2(VP2(:,E(:,1)),VP1(:,E(:,1)))-atan2(VP2(:,E(:,2)),VP1(:,E(:,2))))/2);
  % Even more accurate
  TA = tan( ...
    atan2( ...
      VP1(:,E(:,2)).*VP2(:,E(:,1)) - VP2(:,E(:,2)).*VP1(:,E(:,1)), ...
      VP1(:,E(:,2)).*VP1(:,E(:,1)) + VP2(:,E(:,2)).*VP2(:,E(:,1)))/2);
  W = zeros(size(V,1),size(P,1));
  W(:,E(:,1)) = W(:,E(:,1)) + TA;
  W(:,E(:,2)) = W(:,E(:,2)) + TA;
  W = W./R;
  W = bsxfun(@rdivide,W,sum(W,2));
  if enforce_boundary_conditions
    [b,bc] = boundary_conditions(V,[],P,1:size(P,1),[],E);
    W(b,:) = bc;
  else
    % good luck!
  end
end
