function [VV,Uforward,Ubackward] = untangle(V,F,varargin)
  % UNTANGLE Given a possible self-intersecting mesh (V,F) _untangle_
  % self-intersections by running a (conformalized) mean curvature flow
  % (guaranteed to resolve self-intersections if converges for sphere topology
  % surfaces) and then re-inflating with collision response. This is most
  % similar to the method described in "Consistent Volumetric Discretizations
  % Inside Self-Intersecting Surfaces" [Sacht et al. 2013], which itself is an
  % expansion of the idea in "Interference-Aware Geometric Modeling" [Harmon et
  % al. 2011].
  %
  % VV = untangle(V,F)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   Optional:
  %     'ExpansionEnergy' Followed by one of the following energy types:
  %       'none'  just let eltopo find a valid state and move on
  %       'arap'  gradient descent on ARAP energy after each reverse step.
  %     'FinalEnergy' Followed by one of the following energy types:
  %       'none'  just let eltopo find a valid state and move on
  %       'arap'  gradient descent on ARAP energy after each reverse step.
  %     'Delta' followed by delta value, should roughly be in range
  %       [1e-13,1e13] {1}
  %     'MaxIter' followed by maximum number of iterations {100}
  %     'MaxEnergyIter' followed by maximum number of iterations for energy
  %       minimization {100}
  %     'LaplacianType' followed by 'cotangent' of 'uniform'.
  % Outputs:
  %   VV  #V by dim list of new vertex positions
  %   U  #V by dim by steps list of steps
  %

  function append_Ubackward()
    if main_nargout>=3 
      if isempty(Ubackward)
        Ubackward = zeros([size(VV) 1+(size(Uforward,3)-1)*e_max_eit+f_max_eit]);
      end
      Ubackward(:,:,bit) = VV;
      bit = bit+1;
    end
  end

  delta = 3e-5;
  max_iter = 100;
  laplacian_type = 'cotangent';
  f_max_eit = 100;
  e_max_eit = 10;
  expansion_energy = 'none';
  final_energy = 'none';
  V0 = V;
  rescale_output = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MaxIter','ExpansionEnergy','FinalEnergy','Delta','LaplacianType', ...
      }, ...
    {'max_iter','expansion_energy','final_energy','delta','laplacian_type', ...
      });
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

  if strcmp(expansion_energy,'none')
    e_max_eit = 1;
  end
  if strcmp(final_energy,'none')
    f_max_eit = 1;
  end


  [~,Uforward] = conformalized_mean_curvature_flow( ...
    V,F,'delta',3e-5,'MaxIter',100, ...
    'LaplacianType',laplacian_type, ...
    'UntilSelfIntersectionFree',true, ...
    'RescaleOutput',true);
  
  t = tsurf(F,V,fphong);
  colormap(parula(15));
  axis equal;
  view(2);
  %U = {V};
  %prev_IF = [];
  %while true
  %  first_only = false;
  %  [~,~,IF] = selfintersect(U{end},F,'DetectOnly',true,'FirstOnly',first_only);
  %  title(sprintf('%d',size(IF,1)));
  %  if ~isempty(prev_IF) && size(IF,1) > size(prev_IF,1)
  %    input(sprintf('number of self-intersections is increasing from %d to %d',size(prev_IF,1),size(IF,1)));
  %  end
  %  prev_IF = IF;
  %  if isempty(IF)
  %    break;
  %  end
  %  [U{end+1}] = laplacian_smooth(U{1},F,'uniform',[],0.1,'implicit',U{end},1);
  %  t.Vertices = U{end};
  %  drawnow;
  %end

  Ubackward = [];
  bit = 1;
  main_nargout = nargout;
  VV = Uforward(:,:,end);
  append_Ubackward();

  t.Vertices = VV;
  drawnow;
  for it = size(Uforward,3)-1:-1:0
    if it == 0
      energy  = final_energy;
      max_eit = f_max_eit;
      it = 1;
    else
      energy  = expansion_energy;
      max_eit = e_max_eit;
    end
    switch energy
    case 'arap'
      tol = 1e-7;
      VV_prev = VV;
      delta_t = 0.3;
      for eit = 1:max_eit

        [G,E_prev,R] = arap_gradient(Uforward(:,:,it),F,VV_prev);
        dVV = G;
        while true
          VV = VV_prev - delta_t * dVV;
          [~,E] = arap_gradient(Uforward(:,:,it),F,VV);
          if E<E_prev
            break;
          end
          delta_t = delta_t*0.9;
          assert(delta_t>0);
        end
        [VV,~] = eltopo(VV_prev,F,VV);
        append_Ubackward();
        title(sprintf('step %d/%d, %d/%d',size(Uforward,3)-it,size(Uforward,3)-1,eit,max_eit),'FontSize',15);
        t.Vertices = VV;
        C = normrow(G);
        t.CData = C;
        colorbar;
        caxis([min(C) max(C)]);
        drawnow;
        dVV = max(abs(VV-VV_prev));
        if dVV < tol
          break;
        end
        [~,E] = arap_gradient(Uforward(:,:,it),F,VV);
        if E_prev < E
          delta_t = 0.5*delta_t;
        end
        VV_prev = VV;
      end
    case 'none'
      title(sprintf('step %d/%d',size(Uforward,3)-it,size(Uforward,3)-1),'FontSize',15);
      drawnow;
      [VV,~] = eltopo(VV,F,Uforward(:,:,it));
      append_Ubackward();
      t.Vertices = VV;
    otherwise
      error(sprintf('Unknown energy type (%s)',energy));
    end

    [~,~,IF] = selfintersect(VV,F,'DetectOnly',true,'FirstOnly',true);
    if ~isempty(IF)
      error('eltopo failed');
    end
  end
end
