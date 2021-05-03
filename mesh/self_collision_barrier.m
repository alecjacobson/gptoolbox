function [f,G,H] = self_collision_barrier(V,E,tol)
  % SELF_COLLISION_BARRIER Compute the a barrier function (and its derivatives)
  % for point-edge collisions of a given *strictly feasible** line-complex in
  % 2D. See "Bijective Parameterization with Free Boundaries" [Smith and
  % Schaefer 2015] or for a similar barrier "Incremental Potential Contact:
  % Intersection- and Inversion-free, Large-Deformation Dynamics" [Li et al.
  % 2020].
  % 
  % [f,G,H] = self_collision_barrier(V,E,tol)
  %
  % Inputs:
  %   V  #V by 2 list of input vertex positions
  %   E  #E by 2 list of edge indices into rows of V
  %   tol  distance at which barrier term becomes positive {1e-3}
  % Outputs:
  %   f  scalar total objective value
  %   G  #V*2 gradient 
  %   H  #V*2 by #V*2 sparse Hessian matrix
  %
  % Note: this function uses a special symbolic library trick which generates
  % files self_collision_barrier_cap_sym.m and self_collision_barrier_line_sym.m
  % if they don't already exist. If your change the symbolic-math part of this
  % file, you must delete these automatically generated files.
  %

  function [sqrD,T] = point_segment_squared_distance(P,A,B)
    PA = P-A;
    BA = B-A;
    T = min(max(sum(PA.*BA,2)./sum(BA.^2,2),0),1);
    % vector to closest point
    PC = PA-BA.*T;
    sqrD = sum(PC.^2,2);
  end

  function [f,G,H] = self_collision_barrier_cap(V,E,IJ)
    % vertex i colliding with edge jk

    n = size(V,1);
    % [x x x y y y]
    IJ = [IJ  n+IJ];
    X = V(IJ);

    %% Sanity check
    %p = X(:,1:2:end);
    %c = X(:,2:2:end);
    %pc = p-c;
    %d = sqrt(sum(pc.^2,2));
    %if ~isempty(d) && max(d) > tol && nargout > 1
    %  error
    %end

    % Writing the symbol
    path = mfilename('fullpath');
    aux = [path '_cap_sym.m'];
    % should also check date...
    if ~exist(aux,'file')
      % From here, only touch X
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      sX = sym('X',[1 4]);
      stol = sym('tol',[1 1]);
      sp = sX(1:2:end);
      sc = sX(2:2:end);
      spc = sp-sc;
      sd = sqrt(sum(spc.^2,2));
      sf = barrier(sd,stol);

      hess = @(sf,sX) cell2sym(arrayfun(@(g) gradient(g,sX),gradient(sf,sX),'UniformOutput',false));

      aux_handle = ...
        matlabFunction(sf,gradient(sf,sX),hess(sf,sX),'vars',{sX,stol},'File',aux);
    else
      aux_name = [mfilename('func') '_cap_sym'];
      aux_handle = str2func(aux_name);
    end

    faux=[];gaux = [];Haux = [];
    switch nargout
    case {0,1}
      [faux] = aux_handle(X,tol);
    case 2
      [faux,gaux] = aux_handle(X,tol);
    case 3
      [faux,gaux,Haux] = aux_handle(X,tol);
    end

    % unnecessary indirection?
    f_fun = @(X) faux;
    dfdX_fun = @(X) gaux;
    d2fdX2_fun = @(X) Haux;
    f = sum(f_fun(X));
    if nargout<=1
      return;
    end
    dfdX = dfdX_fun(X);
    G = full_sparse(IJ,ones(size(IJ)),reshape(dfdX,size(IJ)),2*n,1);
    if nargout<=2
      return;
    end
    d2fdX2 = double(d2fdX2_fun(X));


    d2fdX2 = reshape(d2fdX2,[],4*4);
    if psd_project
      d2fdX2 = psd_project_rows(d2fdX2);
    end

    HI = repmat(IJ,[1 1 size(IJ,2)]);
    HJ = permute(repmat(IJ,[1 1 size(IJ,2)]),[1 3 2]);
    H = fast_sparse(HI(:),HJ(:),d2fdX2(:),2*n,2*n);
  end

  function [f,G,H] = self_collision_barrier_line(V,E,IJK)

    n = size(V,1);
    % [x x x y y y]
    IJK = [IJK n+[IJK]];
    X = V(IJK);

    %% Sanity check
    %p = X(:,1:3:end);
    %a = X(:,2:3:end);
    %b = X(:,3:3:end);
    %pa = p-a;
    %ba = b-a;
    %t = sum(pa.*ba,2)./sum(ba.^2,2);
    %pc = pa - ba.*t;
    %d = sqrt(sum(pc.^2,2));
    %if ~isempty(d) && max(d) > tol && nargout > 1
    %  error
    %end

    % Writing the symbol
    path = mfilename('fullpath');
    aux = [path '_line_sym.m'];
    % should also check date...
    if ~exist(aux,'file')
      % From here, only touch X
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      sX = sym('X',[1 6]);
      stol = sym('tol',[1 1]);
      sp = sX(1:3:end);
      sa = sX(2:3:end);
      sb = sX(3:3:end);
      spa = sp-sa;
      sba = sb-sa;
      st = sum(spa.*sba,2)./sum(sba.^2,2);
      spc = spa - sba.*st;
      sd = sqrt(sum(spc.^2,2));
      sf = barrier(sd,stol);

      hess = @(sf,sX) cell2sym(arrayfun(@(g) gradient(g,sX),gradient(sf,sX),'UniformOutput',false));
      aux_handle = ...
        matlabFunction(sf,gradient(sf,sX),hess(sf,sX),'vars',{sX,stol},'File',aux);
    else
      aux_name = [mfilename('func') '_line_sym'];
      aux_handle = str2func(aux_name);
    end

    faux=[];gaux = [];Haux = [];
    switch nargout
    case {0,1}
      [faux] = aux_handle(X,tol);
    case 2
      [faux,gaux] = aux_handle(X,tol);
    case 3
      [faux,gaux,Haux] = aux_handle(X,tol);
    end

    % unnecessary indirection?
    f_fun = @(X) faux;
    dfdX_fun = @(X) gaux;
    d2fdX2_fun = @(X) Haux;
    f = sum(f_fun(X));
    if nargout<=1
      return;
    end
    dfdX = dfdX_fun(X);
    G = full_sparse(IJK,ones(size(IJK)),reshape(dfdX,size(IJK)),2*n,1);
    if nargout<=2
      return;
    end
    d2fdX2 = double(d2fdX2_fun(X));

    d2fdX2 = reshape(d2fdX2,[],6*6);
    if psd_project
      d2fdX2 = psd_project_rows(d2fdX2);
    end

    HI = repmat(IJK,[1 1 size(IJK,2)]);
    HJ = permute(repmat(IJK,[1 1 size(IJK,2)]),[1 3 2]);
    H = fast_sparse(HI(:),HJ(:),d2fdX2(:),2*n,2*n);

  end

  psd_project = true;
  % The max is not needed because we're handling that explicitly
  % 
  % Smith and Schaefer
  %barrier = @(d,tol) (tol./d - 1).^2;
  % IPC
  barrier = @(D,tol) -(D-tol).^2.*log(D./tol);

  b = unique(E);
  [b1,b2] = box_each_element(V,b);
  [E1,E2] = box_each_element(V,E);
  I = box_intersect(b1-tol,b2+tol,E1-tol,E2+tol);
  I = I(~any(b(I(:,1))==E(I(:,2),:),2),:);
  BC = barycenter(V,E);

  [sqrD,T] = point_segment_squared_distance(V(b(I(:,1)),:),V(E(I(:,2),1),:),V(E(I(:,2),2),:));

  keep = find(sqrD<tol^2);

  I = I(keep,:);
  T = T(keep);

  %BC = barycenter(V,E);
  %plot_edges(V,E,'-ok');
  %hold on;
  %sct(V(b(I(:,1)),:),'or','LineWidth',2);
  %I(:,2)
  %size(E)
  %sct(BC(I(:,2),:),'og','LineWidth',2);
  %hold off;

  cap = T==0 | T==1;
  IJKline = [b(I(~cap,1)) E(I(~cap,2),:)];

  sqrD = sqrD(keep);

  Icap = I(cap,:);
  Tcap = T(cap);
  Jcap = E(Icap(:,2),1);
  Jcap(Tcap==1) = E(Icap(Tcap==1,2),2);
  IJcap = unique(sort([b(Icap(:,1)) Jcap],2),'rows');

  switch nargout
  case {0,1}
    [fline] = self_collision_barrier_line(V,E,IJKline);
    [fcap] = self_collision_barrier_cap(V,E,IJcap);
    f = fline + fcap;
  case 2
    [fline,Gline] = self_collision_barrier_line(V,E,IJKline);
    [fcap,Gcap] = self_collision_barrier_cap(V,E,IJcap);
    f = fline + fcap;
    G = Gline + Gcap;
  case 3
    [fline,Gline,Hline] = self_collision_barrier_line(V,E,IJKline);
    [fcap,Gcap,Hcap] = self_collision_barrier_cap(V,E,IJcap);
    f = fline + fcap;
    G = Gline + Gcap;
    H = Hline + Hcap;
  end



end
