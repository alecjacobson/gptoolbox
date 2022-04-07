function [UU,rU,rF] = distribute_on_ground(VV,FF,varargin)
  % [UU,rU,rF] = distribute_on_ground(VV,FF,'ParameterName',ParameterValue, ...)
  %
  % Given a list of meshes. Iteratively rest each mesh on the ground, ensuring
  % that new meshes don't intersect with those already placed.
  % 
  % Inputs:
  %   VV  #VV-long list of mesh vertex positions arrays
  %   FF  #VV-long list of mesh triangle index arrays
  %   Optional:
  %     'Fast'  followed by whether to use fast overly conservative bounding box
  %     intersection tests instead of mesh-based ones.
  % Outputs:
  %   UU  #VV-long list of new mesh vertex positions
  %   rU  #rU-long concatenation of all new meshes
  %   rF  #rF-long triangle mesh indices into rows of rU
  %
  % Known issues: this has only been tested for models that fit more or less in
  % the unit sphere). 
  %

  % Map of parameter names to variable names
  fast = false;
  params_to_variables = containers.Map( ...
    {'Fast'}, ...
    {'fast'});
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

  vol = zeros(numel(FF),1);
  for i = 1:numel(FF)
    [~,vol(i)] = centroid(VV{i},FF{i});
  end
  [~,sI] = sort(vol,'descend');
  VV = VV(sI);
  FF = FF(sI);

  UU = VV;
  rU = [];
  rF = [];
  for i = 1:numel(FF)
    UU{i} = rest_on_ground(UU{i},FF{i});
    best_Ui = UU{i};
    best_d = inf;
    if i > 1
      if fast
        rB1 = min(rU);
        rB2 = max(rU);
      end
      while isinf(best_d)
        tries = 10 + fast*100;
        for t = 1:tries
          cen = [randn(1,2) 0];
          R2 = rand_rotation(2);
          R = blkdiag(R2,1);
          R = R*det(R);
          R = eye(3,3);
          Ui = UU{i}*R+cen;
          d = sum((Ui-UU{i}).^2,2);
          if fast
            B1 = min(Ui);
            B2 = max(Ui);
            IF = box_intersect(rB1,rB2,B1,B2);
          else
            IF = intersect_other(rU,rF,Ui,FF{i},'FirstOnly',true);
          end
          if isempty(IF)
            if d < best_d
              best_Ui = Ui;
              best_d = d;
            end
          end
        end
      end
    end
    if ~fast
      rF = [rF;size(rU,1)+FF{i}];
    end
    rU = [rU;best_Ui];
    UU{i} = best_Ui;
  end
  % Unsort output
  UU(sI) = UU;
end
