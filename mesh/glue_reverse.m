function [VV,FF,RIM] = glue_reverse(V,F,varargin)
  % GLUE_REVERSE Glue a mesh to itself at its boundary reversing the
  % z-coordinates of one of the copies. Useful for creating closed meshes from
  % inflated disks.
  %
  % [VV,FF,RIM] = glue_reverse(V,F)
  %
  % Inputs:
  %   V  #V by 3  list of vertex positions
  %   F  #F by 3  list of face indices
  %   Optional:
  %     'Boundary'  followed #E by 2 by subset of the boundary edges to glue
  %       along. {all boundary edges}
  % Outputs:
  %  VV #VV by 3 list of new vertex positions
  %  FF #FF by 3 list of new face indices
  %  RIM #VV by 1 list of reindexing indices such that VV = V(RIM,:)
  %

  % default values
  out = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Boundary'}, ...
    {'out'});
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

  % number of domain vertices
  n = size(V,1);
  dim = size(V,2);
  % concatenate vertex lists with Z-coord of copy flipped
  VV = [V; V(:,1:2) -1*V(:,3:dim)];
  % get list of unique boundary vertex indices
  if isempty(out)
    out = unique(reshape(outline(F),[],1));
  end
  VV(out,3:dim) = 0;
  % get list of new faces as if we'd just concatenate vertex lists, normals of new
  % faces are flipped
  Fflip = fliplr(F);
  % Tell new boundary vertices to point to original ones
  Fflip(~ismember(Fflip,out)) = n+Fflip(~ismember(Fflip,out));
  FF = [F;Fflip];
  % rearrange vertices so those used by faces come first
  [VV,~,FF,IM] = faces_first(VV,[],FF);
  % only keep up to those used by faces (throw out boundary copies)
  VV = VV(1:max(FF(:)),:);
  %IM = IM(1:max(FF(:)),:);
  RIM = zeros(size(VV,1),1);
  RIM(IM) = [1:n 1:n];
  RIM = RIM(1:max(FF(:)),:);
end
  
