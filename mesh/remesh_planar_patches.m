function [W,G,IM,C] = remesh_planar_patches(V,F,varargin)
  % REMESH_PLANAR_PATCHES  Find nearly planar patches and retriangulate them.
  %   (V,F) should probably not be self-intersecting, at least not near the
  %   planar patches. This will attempt to maintain non-manifold edges
  %   (untested).
  %
  % [W,G] = remesh_planar_patches(V,F)
  % [W,G,IM,C] = remesh_planar_patches(V,F,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 3 list of vertices
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'Force' followed by whether to remesh even if it will not improve
  %        the number of vertices {false}.
  %     'MinSize' followed by minimum number of facets in group to consider for
  %       remeshing {4}
  %     'MinDeltaAngle' followd by minimum change in angle between neighboring
  %     facets to be considered co-planar {pi-1e-5}
  %     'Except' followed by list of faces not to remesh
  %     'TriangleFlags' followed by string of flags to pass to triangle
  %
  % Outputs:
  %   W  #W by 3 list of output mesh positions
  %   G  #G by 3 list of output triangle indices
  %   IM  indices from final remove_unreferenced
  %   C  #F list of patch indices
  %

  % force remesh even if there are no internal vertices to remove
  force_remesh = false;
  % Minmimum angle for two neighboring facets to be considered coplanar
  % in radians
  min_delta_angle = pi-1e-5;
  %min_delta_angle = pi-1e-3;
  min_size = 4;
  triangle_flags = '';
  except = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Force','MinSize','MinDeltaAngle','TriangleFlags','Except'}, ...
    {'force_remesh','min_size','min_delta_angle','triangle_flags','except'});
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

  % The right way to deal with them would be to make a list of feature edges
  % (aka "exterior edges") on the original mesh and be sure these end up in any
  % planar pataches. They would include the patch boundaries and any other
  % internal non-manifold edges.
  NME = nonmanifold_edges(F);
  C = planar_patches(V,F,'MinDeltaAngle',min_delta_angle,'Except',except);

  [UC,~,ic] = unique(C);
  ucounts = histc(C,UC);
  counts = ucounts(ic);

  %tsurf(F,V,'CData',C);
  %axis equal;
  %drawnow;

  W = V;
  G = F(counts<=min_size,:);
  % loop over connected components
  %for c = reshape(UC(ucounts>min_size),1,[])
  for c = 1:max(C)
    if ucounts(c) > min_size
      Fc = F(C==c,:);
      assert(size(Fc,1)>min_size);
      %S = statistics(V,Fc,'Fast',true);
      %% Not dealing with holes yet
      %if S.num_boundary_loops > 1
      %  warning('Skipping high-genus planar patch...');
      %  G = [G;Fc];
      %  continue;
      %end
      %assert(S.num_boundary_loops == 1);
      uf = unique(Fc);
      % Fit plane to all points
      [~,A] = affine_fit(V(uf,:));
      % Only keep outline points
      O = outline(Fc);

      % Non manifold edges on this patch
      NMEc = NME(all(ismember(NME,uf),2),:);
      E = [O;NMEc];

      uo = unique(E);
      % Nothing to gain by remeshing
      if ~force_remesh && numel(uo) == numel(uf)
        G = [G;Fc];
        continue;
      end
      % only boundary edges projected to plane
      Vuo = V(uo,:)*A;
      J = 1:size(V,1);
      J(uo) = 1:numel(uo);
      % Remap E to Vu
      O = J(O);
      E = J(E);

      % THIS IS INSANELY SLOW
      % DT = delaunayTriangulation(Vuo(:,1),Vuo(:,2));
      % DT.Constraints = E;
      % Gc = DT.ConnectivityList;
      % Triangle is way faster...
      %[Wuo,Gc] = triangle(Vuo,E,[],'Quiet','Flags',triangle_flags);
      [Wuo,Gc] = triangulate(Vuo,E,'Flags',triangle_flags);

      %tsurf(Gc,Wuo);
      %input('');
      assert(size(Gc,1) >= 1);
      assert(size(Wuo,1) >= size(Vuo,1));

      % easier to remove holes post hoc than pass hole positions to triangle
      % However, we should call triangulate again with the holes removed...
      AE = adjacency_matrix(E);
      [~,CE] = conncomp(AE);
      if max(CE) > 1
        w = winding_number(Vuo,O,barycenter(Wuo,Gc));
        % should only be 1s and 0s
        Gc = Gc(abs(w)>0.1,:);
      end


      % remap Gc to V
      uo = [uo;(size(W,1)+(1:(size(Wuo,1)-size(Vuo,1))))'];
      Wnew = Wuo((size(Vuo,1)+1):end,:);
      VA = V*A;
      %tsurf(Fc,VA);
      %hold on;
      %scatter(Wnew(:,1),Wnew(:,2));
      %hold off;
      %pause

      I = in_element(VA,Fc,Wnew);
      [~,I] = find(I);
      assert(numel(I) == size(Wnew,1));
      B = barycentric_coordinates( ...
        Wnew,VA((Fc(I,1)),:),VA((Fc(I,2)),:),VA((Fc(I,3)),:));
      Wnew = sum(bsxfun(@times, permute(B,[1 3 2]), ...
        cat(3,V(Fc(I,1),:),V(Fc(I,2),:),V(Fc(I,3),:))),3);
      W = [W;Wnew];
      Gc = uo(Gc);
      
      %tsurf(Gc,W);
      %pause
      % old mean normal
      oldN = mean(normalizerow(normals(V,Fc)));
      N = mean(normalizerow(normals(W,Gc)));
      flip  = dot(oldN,N) < 0;
      if flip
        Gc = fliplr(Gc);
      end
      %tsurf(Gc,W);
      %input('');
      G = [G;Gc];
      assert(max(G(:))<=size(W,1));
    end
  end

  [W,IM] = remove_unreferenced(W,G);
  G = IM(G);

end
