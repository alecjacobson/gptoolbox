function D = is_extrinsic_delaunay(V,F,varargin)
  % IS_EXTRINSIC_DELAUNAY Determine if each edge in a 2D mesh (V,F) is *extrinsically*
  % Delaunay. 
  %
  % D = is_extrinsic_delaunay(V,F)
  % D = is_extrinsic_delaunay(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangles indices
  %   Optional:
  %     'Tol'  followed by tolerance circumradius check, in length units
  % Outputs:
  %   D  #F by 3 list of bools revealing whether opposite edges are locally
  %     Delaunay. Boundary edges are by definition Delaunay.
  %

  % default values
  tol = 0;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol'}, {'tol'});
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

  assert(size(V,2) == 2,'V should be 2D');
  O = outline(F);
  [~,B] = on_boundary(F);

  % append point at infinity to deal with Meshes with boundaries
  if ~isempty(O)
    FF = [F;repmat(size(V,1)+1,size(O,1),1) O(:,[2 1])];
    VV = [V;inf inf];
    DD = is_extrinsic_delaunay(VV,FF,varargin{:});
    assert(all(size(DD)==size(FF)));
    D = DD(1:size(F,1),:);
    D(B) = true;
    return;
  end

  % Mesh better be manifold
  [Fp, Fi] = triangle_triangle_adjacency(F);
  % Mesh better not have a boundary (handled above)
  assert(all(Fp(:)>0));
  assert(all(Fi(:)>0));
  opp_ind = sub2ind(size(F),Fp,Fi);
  Fopp = F(opp_ind);
  [R,C] = circumradius(V,F);
  D = reshape(normrow(V(Fopp(:),:)-repmat(C,3,1))-repmat(R,3,1)>=-abs(tol),size(F));
  % symmetrize 
  D = D & D(opp_ind);

end
