function [IO,A] = apply_texture_map(tsh,UV,T,varargin);
  % [IO,A] = apply_texture_map(tsh,UV,T)
  % 
  % Inputs:
  %   tsh  handle to trisurf
  %   UV  #tsh.Vertices by 2 list of UV-coordinates into [1:size(T,2) by 1:size(T,1)]
  %   T  tc-height by tc-width by 3 rgb texture image
  % Output:
  %   IO  gcf-height by gcf-width by 3 rgb image of current frame with mesh
  %     but face colors are replaced by texture color
  %   A  gcf-height by gcf-widsth boolean image of where the mesh (faces) are.
  %
  % Known Issues:
  %   UVs here are mapping to _rectangle_ of TC rather than the usual unit
  %   square. This is to facilitate use with bwmesh.  If reading .obj + .png,
  %   this should be adapted and possibly controled with a flag.
  %
  interp_args = {};
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'InterpArgs'}, ...
    {'interp_args'});
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

  [I,A,FI,B] = shader(tsh);
  [II,~,IV] = find(FI(:));
  B1 = B(:,:,1);B2 = B(:,:,2);B3 = B(:,:,3);
  F = tsh.Faces;
  PUV = zeros(numel(B1),2);
  PUV(II,:) =  ...
    UV(F(IV,1),:).*B1(II) + UV(F(IV,2),:).*B2(II) + UV(F(IV,3),:).*B3(II);
  PUV = reshape(PUV,[size(B,1),size(B,2),size(PUV,2)]).*A;
  [XT,YT] = meshgrid(1:size(T,2),size(T,1):-1:1);
  IT = [];
  T = im2double(T);
  IT(:,:,1) = interp2(XT,YT,T(:,:,1),PUV(:,:,1),PUV(:,:,2),interp_args{:}).*A;
  IT(:,:,2) = interp2(XT,YT,T(:,:,2),PUV(:,:,1),PUV(:,:,2),interp_args{:}).*A;
  IT(:,:,3) = interp2(XT,YT,T(:,:,3),PUV(:,:,1),PUV(:,:,2),interp_args{:}).*A;
  IT(isnan(IT(:))) = 0;
  IO = IT.*A + I.*~A;
end
