function [A,R,J,I,W,G1,S,C] = spline_G1_subspace(P,C,varargin)
  % [A,R,J,I,W,G1,S] = spline_G1_subspace(P,C)
  % 
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  % Optional:
  %   'Tol' followed by the tolerance used to determine sharp corners in
  %       input in radians [0,π] {0.0447}
  % Outputs:
  %   A  #P by #R sparse matrix mapping reduced coordinates (R) to full, input
  %     coordinates so that all C0 and G1 constraints are satisfied
  %   R  #R by dim list of reduced control point locations so that A*R ≈ P
  %     (approximate because of `tol`).
  %   J  #R list of indices so that R = P(J,:)
  %   I  #I list of indices into rows of P where C0 or G1 constraints were enforced
  %   W  #I by 2 list of pairs of indices into rows of  C where C0 or G1
  %     constraints were enforced.
  %   G1  #G1 list of indices into I where G1 enforced
  %   S  #G1 list of scale factors indicating S T_W(G1,2) = T_W(G1,1)
  %   C  #C by 4 possibly flipped versions of C to match tangent directions used
  %     to compute S
  %   
  % Example:
  %   [P,C] = cubic_subdivide([0 0;1 1;2 -1;3 0],[0.25 0.5 0.75]);
  %   A = spline_G1_subspace(P,C);
  %   % Random perturbation
  %   Pr = P + 0.1*randn(size(P));
  %   clf;
  %   hold on;
  %   plot_spline(Pr,C,'Color',blue);
  %   % Project onto space of G1 continuous curves
  %   plot_spline(A*((A'*A)\(A'*Pr)),C,'Color',orange);
  %   hold off;
  %

  tol = 0.0447;
  W = [];
  G1 = [];
  I = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol','W','G1','I'}, ...
    {'tol','W','G1','I'});
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
  assert(tol >=0 & tol <= pi);

  % The following will be destructive to C
  OC = C;

  if isempty(W) 
    E = [C(:,1) C(:,end)];
    [K,A] = manifold_patches(E);
    % O(#components)
    total_num_fit_subsequences = 0;
    I = [];
    W = [];
    for k = 1:max(K)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % factor out this O(#comps #curves) by pre-sorting etc.
      keep = find(K==k);
    
      if numel(keep)>1
        Ek = [C(keep,1) C(keep,end)];
        if numel(keep) == 2
          % special handling for loop with two segments since this is degenerate
          % when treated as an undirected graph
          if Ek(2,1) == Ek(1,2)
            Ik = [Ek(1,1);Ek(1,2);Ek(2,2)];
            F = [1;1];
          elseif Ek(2,2) == Ek(1,1)
            Ik = [Ek(1,2);Ek(1,1);Ek(2,1)];
            F = [2;2];
          else
            assert(Ek(2,2) == Ek(1,2));
            Ik = [Ek(1,1);Ek(1,2);Ek(2,1)];
            F = [1;2];           
          end
          J = [1;2];
        else
          [Ik,J,F] = edges_to_path(Ek);
        end
        C(keep(F==2),:) = fliplr(C(keep(F==2),:));
        is_loop = Ik(1) == Ik(end);
        if is_loop
          Wk = [J(1:end) J([2:end 1])];
          Ik = Ik(2:end);
        else
          Wk = [J(1:end-1) J([2:end])];
          Ik = Ik(2:end-1);
        end
        assert(size(Wk,1) == numel(Ik));
        assert(all(C(keep(Wk(:,1)),end) == C(keep(Wk(:,2)),1)));
        I = [I;Ik];
        W = [W keep(Wk)'];
      end
    
    end
    W = W';
  end

  if isempty(I)
    I = C(W(:,1),end);
  end

  assert(all(C((W(:,1)),end) == C((W(:,2)),1)));
  assert(all(C((W(:,1)),end) == I));
  
  Tin =  P(C(W(:,1),end),:) - P(C(W(:,1),end-1),:);
  Tout = P(C(W(:,2),2),:) - P(C(W(:,2),1),:);
  if isempty(G1)
    % Every row in W is a candidate for a G1 constraint
    G1 = find(acos(sum(normalizerow(Tin).*normalizerow(Tout),2)) < tol);
  elseif isscalar(G1)
    G1 = repmat(G1,size(W,1),1);
  end
  S = normrow(Tout(G1,:))./normrow(Tin(G1,:));
  
  % Might be faster to build transpose: `any` below is fast though
  A = speye(size(P,1),size(P,1));
  % C₂² = C₄¹ - S (C₃¹ - C₄¹)
  % C₂² = (1+S) C₄¹ - S C₃¹
  % C₂² = -S C₃¹ + (1+S) C₄¹
  A( C(W(:,2),1), :) = sparse(1:size(W,1),C(W(:,1),end),1,size(W,1),size(P,1));
  A( C(W(G1,2),2), :) = sparse( ...
    repmat(1:numel(S),2,1)', ...
    C(W(G1,1),[end-1 end]), ...
    [-S (1+S)], ...
    numel(S),size(P,1));
  J = find(any(A));
  A = A(:,J);
  R = P(J,:);

end
