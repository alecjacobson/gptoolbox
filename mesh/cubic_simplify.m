function [Pf,Cf,total_num_fit_subsequences] = cubic_simplify(P,C,varargin)
  % [Pf,Cf,total_num_fit_subsequences] = cubic_simplify(P,C,varargin)
  % 
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  %   Optional:
  %     'FlatTol' followed by the tolerance for cubic_flat_eval {0.1}
  %     'FitTol'  followed by the tolerance of the fit {1}
  %     'MaxLength'  followed by the maximum length between samples {1}
  %     'G1Tol'  followed by the tolerance used to determine sharp corners in
  %       input in radians [0,π] {0.0447}
  % Outputs:
  %   Pf  #Pf by dim list of control point locations
  %   Cf  #Cf by 4 list of indices into Pf of cubic Bezier curves
  % 
  flat_tol = 0.1;
  fit_tol = 1e0;
  max_length = 1;
  G1_tol = 0.0447;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'FlatTol','FitTol','G1Tol','MaxLength'}, ...
    {'flat_tol','fit_tol','G1_tol','max_length'});
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

  assert(G1_tol >=0 & G1_tol <= pi);

  bbd = normrow(max(P)-min(P));


  [P,~,~,C] = remove_duplicate_vertices(P,0,'F',C);
  E = [C(:,1) C(:,end)];
  [K,A] = manifold_patches(E);
  % ignore orientation and connect into strips
  Pf = [];
  Cf = [];
  CC = {};
  % O(#components)
  total_num_fit_subsequences = 0;
  for k = 1:max(K)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % factor out this O(#comps #curves) by pre-sorting etc.
    keep = K==k;
    Ek = E(keep,:);
    Ck = C(keep,:);


    % find path along segements and reorient if necessary
    [I,J,F] = edges_to_path(Ek);
    Ck = Ck(J,:);
    % Reorient
    Ck(F==2,:) = fliplr(Ck(F==2,:));
    is_loop = I(1) == I(end);
    if is_loop
      W = [1:size(Ck,1);2:size(Ck,1) 1]';
    else
      W = [1:size(Ck,1)-1;2:size(Ck,1)]';
    end
    assert(all(Ck(W(:,1),end)==Ck(W(:,2),1)));
    Tin =  normalizerow(P(Ck(W(:,1),end),:) - P(Ck(W(:,1),end-1),:));
    Tout = normalizerow(P(Ck(W(:,2),2),:) - P(Ck(W(:,2),1),:));
    G1 = acos(sum(Tin.*Tout,2)) < G1_tol;
    if any(G1) && ~all(G1) && is_loop
      % rotate so that a segment ending in a non-G1 point is last
      i = find(~G1,1,'last');
      Ck = Ck([i+1:end,1:i],:);
      G1 = G1([i+1:end,1:i]);
    end
    CS = [0;find(~G1)];
    if CS(end) ~= size(Ck,1)
      CS(end+1) = size(Ck,1);
    end
    % For each piecewise G1 subsequence
    Pkf = [];
    Ckf = [];
    total_num_fit_subsequences = total_num_fit_subsequences + numel(CS)-1;
    for s = 1:numel(CS)-1
      if CS(s+1)-CS(s)>1
        % Don't use this because it's being clever about duplicate vertices
        %[V,E,I] = spline_to_poly(P,Ck(CS(s)+1:CS(s+1),:),flat_tol);
        V = [];
        for c = CS(s)+1:CS(s+1);
          [Vc,Tc] = cubic_flat_eval(P(Ck(c,:),:),flat_tol);
          VTc = [Vc Tc];
          Ec = [1:size(Vc,1)-1;2:size(Vc,1)]';
          [VTc,Ec] = upsample(VTc,Ec,'Iterations',inf,'OnlySelected',@(V,E) edge_lengths(V(:,1:end-1),E)>max_length);
          Tc = VTc(:,end);
          Vc = cubic_eval(P(Ck(c,:),:),Tc);
          I = edges_to_path(Ec);
          Vc = Vc(I,:);
  
          if c ~= CS(s+1)
            Vc = Vc(1:end-1,:);
          end
          V = [V Vc'];
        end
        V = V';
        perturbed_last = false;
        if is_loop && numel(CS)-1 == 1 && any(G1)
          perturbed_last = true;
          V(end,:) = V(end,:) + bbd*1e-10;
        end
        V = V([1;any(diff(V),2)]~=0,:);
        Ps = cell2mat(fit_cubic_bezier(V,fit_tol));
        if perturbed_last
          Ps(end,:) = Ps(1,:);
        end
        [Ps,~,Cs] = remove_duplicate_vertices(Ps,0);
        Cs = reshape(Cs,4,[])';

        %clf;
        %plt(V(:,1:2));
        %hold on;
        %[V,Ew] = spline_to_poly(P(:,1:2),Ck,1);
        %tsurf(Ew,V,'EdgeColor','r');
        %plot_spline(Ps(:,1:2),Cs);
        %hold off;
        %axis equal;
        %pause

      end
      if CS(s+1)-CS(s) == 1 || size(Cs,1) >= CS(s+1)-CS(s)
        % we shit the bed and made things worse.
        Cs = Ck(CS(s)+1:CS(s+1),:);
        [Us,~,Cs] = unique(Cs(:));
        Cs = reshape(Cs,[],4);
        Ps = P(Us,:);
      end
      Ckf = [Ckf size(Pkf,2)+Cs'];
      Pkf = [Pkf Ps'];
    end
    Pkf = Pkf';
    Ckf = Ckf';
    [Pkf,~,~,Ckf] = remove_duplicate_vertices(Pkf,0,'F',Ckf);
    Cf = [Cf size(Pf,2) + Ckf'];
    Pf = [Pf Pkf'];
  end
  Cf = Cf';
  Pf = Pf';
  [Pf,~,~,Cf] = remove_duplicate_vertices(Pf,0,'F',Cf);

end
