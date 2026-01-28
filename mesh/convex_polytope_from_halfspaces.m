function [pV,pPI,pPC,dV,dF,dN,x0] = convex_polytope_from_halfspaces(P,x0_in)
  % [pV,pPI,pPC,dV,dF,dN,x0] = convex_polytope_from_halfspaces(P,x0)
  %
  % Construct the convex polytope formed by intersecting halfspaces of the form:
  %
  % pᵢ⋅[x 1] ≤ 0
  % 
  % Inputs:
  %   P  #P by 4 list of halfspace equation coefficients
  %   x0  1 by ndim list of interior point (optional)
  % Outputs:
  %   pV  #pV by ndim list of primal points (intersection of halfspaces)
  %   pPI  #pPI stream of polygon indices into rows of pV
  %   pPC  #pPC list of cumulative sum of polygon valences
  %   dV  #dV by ndim list of dual points 
  %   dF  #dF by 3 list of triangle that triangulate dual faces
  %   dN  #dN by ndim list of normals to the dual faces (in dual space)
  %   x0  1 by ndim list of interior point (if not provided, computed)
  % 
  % Example:
  %   [pV,pPI,pPC,dV,dF,dN] = convex_polytope_from_halfspaces(P,x0);
  %   % Better to snap in dual space (less arithmetic has been done)
  %   %[sP,uI,uJ,sF] = remove_duplicate_vertices(pV,2^-40,'F',pF);
  %   [~,I,J] = unique(round(dN*2^40), 'rows');
  %   [pF,pI,pE] = polygons_to_triangles(pPI,pPC);
  %   sP = pV(I,:);
  %   sF = J(pF);
  %   keep = sF(:,1)~=sF(:,2) & sF(:,2)~=sF(:,3) & sF(:,3)~=sF(:,1);
  %   sF = sF(keep,:);
  %   sI = pI(keep);
  % 

  ndim = size(P,2)-1;
  n_planes = size(P,1);
  if n_planes < ndim+1
      error('convex_polytope_from_halfspaces:Infeasible', ...
            'No solution found. The halfspaces may not form a convex polytope.');
  end
  a = P(:, 1:ndim);
  b = P(:, ndim+1);
  % Clean this up with compute_x0_as_chebyshev_center below
  x0 = [];
  if nargin>=2 && ~isempty(x0_in)
    if max(P*[x0_in 1]') < -1e-13
      % accept the input
      x0 = x0_in;
    end
  end

  compute_x0_as_chebyshev_center = isempty(x0);
  %% Should this be switched to Chebeshev center? in text_x0?
  %c = zeros(ndim+1,1);
  %c(end) = -1;
  %a_ub = [a ones(n_planes,1)];
  %b_ub = -b;
  %lb = -inf(ndim+1, 1);
  %lb(end) = 0;
  %ub = inf(ndim+1, 1);
  %%sol = linprog(c, a_ub, b_ub, [], [], lb, ub);
  %% Quiet
  %options = optimoptions('linprog', 'Display', 'off');
  %sol = linprog(c, a_ub, b_ub, [], [], lb, ub, options);
  %if isempty(sol)
  %  % error('No solution found. The halfspaces may not form a convex polytope.');
  %  % Throw error with code
  %  error('convex_polytope_from_halfspaces:Infeasible', ...
  %        'No solution found. The halfspaces may not form a convex polytope.');
  %end
  %x0 = reshape(sol(1:ndim), 1,ndim);

  success = false;
  method = 'convhulln';
  while ~success
    if compute_x0_as_chebyshev_center
      x0 = chebyshev_center(P);
      if isempty(x0)
        % error('No solution found. The halfspaces may not form a convex polytope.');
        % Throw error with code
        %
        % I claim if cond(P) or cond(P(:,1:ndim)) is large then the enclosed
        % volume is near zero.
        %
        % Same issue below. How to update the connectivity? Merge removed face
        % with closest plane?
        error('convex_polytope_from_halfspaces:Infeasible', ...
              'Unable to find interior point.');
      end
    end
    dV = -a./(sum(a.*x0,2)+b);
    if any(isnan(dV(:))) || any(isinf(dV(:)))
      % Q: Is it safe to assume the presence of an Inf in dV means the dual volume is
      % infinite and therefore the primal volume is zero?
      % H: Maybe it's sufficient but doesn't seem necessary.
      %
      % Even so there's still the question of defining the updated connectivity.
      %[x0,lp] = chebyshev_center(P);
      error('convex_polytope_from_halfspaces:Infeasible', ...
            'Degenerate dual vertex positions.');
    end
    % This always gives triangle faces even if there are polygons. I claim this is
    % fine because we'll deal with merging later (see example).
    %
    % Careful, convhulln flips orientation
    switch method
    case 'convhull'
      try
        dF = convhull(dV);
        success = true;
      catch ME
        if strcmp(ME.identifier,'MATLAB:convhull:NotEnoughPtsConvhullErrId')
          if compute_x0_as_chebyshev_center
            % retry with convhulln
            method = 'convhulln';
          else
            % sometimes convhull(dV) fails and using a better x0 would work.
            warning('had to compute x0');
            compute_x0_as_chebyshev_center = true;
          end
        end
      end
    case 'convhulln'
      % Last resort.
      warning('had to use convhulln');
      dF = fliplr(convhulln(dV));
      success = true;
    end
  end

  switch ndim
  case 2
    dF = [dF(1:end-1) dF([2:end-1 1])];
    starts = zeros(size(dF,1),1);
    starts(dF(:,1)) = 1:size(dF,1);
    pE = [(1:size(dF,1))' starts(dF(:,2))];
    % idiotic unrolling
    pPI = reshape(pE',[],1);
    pPC = cumsum([0;repmat(2,size(pE,1),1)]);
  case 3
    [~,pPI,pPC] = dual(dV,dF);
  end

  % Don't snap primal points just yet
  I = (1:size(dF,1))';
  % barycenter to find point on each plane
  BC = barycenter(dV,dF(I,:));
  dx0 = zeros(1,size(BC,2));
  % offset of plane to dx0
  dN = normalizerow(normals(dV,dF));
  c = sum(dN(I,:).*BC,2);
  pV = x0 + dN(I,:)./c;

end
