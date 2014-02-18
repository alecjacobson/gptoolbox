function FF = flip_ears(V,F,varargin)
  % FLIP_EARS Flip ears (trianles with two boundary edges in a mesh 
  %
  % FF = flip_ears(V,F)
  % FF = flip_ears(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertices
  %   F  #F by 3 list of triangle mesh indices
  %   Optional:
  %     'PlanarOnly' followed by whether to flip only if ears and their flops
  %       form a planar quad. {false}
  %     'PlanarEpsilon' Epsilon used to determine planarity {1e-8}.
  % Outputs:
  %   FF  #F by 3 list of new triangle mesh indices
  %



  % default values
  planar_only = false;
  epsilon = 1e-8;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PlanarOnly','PlanarEpsilon'}, ...
    {'planar_only','epsilon'});
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
  
  % Must be manifold (actually only ears need to be)
  [Fp, Fi] = tt(F);
  ears = find(sum(Fp==-1,2)==2);
  [flops,ear_opp] = max(Fp(ears,:),[],2);
  flop_opp = Fi(sub2ind(size(Fp),ears,ear_opp));
  % should check that ears are not neighbors
  
  % Check that ears + flops are planar 
  
  % Quad
  %Q = [F(ears,:) F(sub2ind(size(F),flops,flop_opp))];
  if planar_only
    V(:,end+1:3) = 0;
    Near = normalizerow(normals(V,F(ears,:)));
    Nflop = normalizerow(normals(V,F(flops,:)));
    D = (1-sum(Near.*Nflop,2))<epsilon;
    % only keep planar pairs
    ears = ears(D);
    flops = flops(D);
    ear_opp = ear_opp(D);
    flop_opp = flop_opp(D);
  end
  
  
  % output
  FF = F;
  % ear_opp+2 --> flop_opp
  FF(sub2ind(size(F),ears,mod(ear_opp+1,3)+1)) = F(sub2ind(size(F),flops,flop_opp));
  FF(sub2ind(size(F),flops,mod(flop_opp+1,3)+1)) = F(sub2ind(size(F),ears,ear_opp));
  
end
