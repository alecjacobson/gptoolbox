function FF = flip_ears(V,F,varargin)
  % FLIP_EARS Flip ears (triangles with two boundary edges)
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
  %     'FlipAndClip' followed by whether to Flip ears then clip any new ears,
  %       then repeat until convergence {false}
  % Outputs:
  %   FF  #F by 3 list of new triangle mesh indices
  %
  % Known bugs: (V,F) mmust be manifold near ears.
  %



  % default values
  planar_only = false;
  epsilon = 1e-8;
  flip_and_clip = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PlanarOnly','PlanarEpsilon','FlipAndClip'}, ...
    {'planar_only','epsilon','flip_and_clip'});
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


  % output
  FF = F;
  while true
    F = FF;
    [ears,ear_opp,flops,flop_opp] = find_ears(F);
    % should check that ears are not neighbors, and then do what?

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

    if flip_and_clip
      [ears,ear_opp,flops,flop_opp] = find_ears(FF);
      if isempty(ears)
        break;
      else
        FF(ears,:) = [];
      end
    else
      break;
    end
  end

end
