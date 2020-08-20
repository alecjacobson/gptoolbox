function [VV,FF,J,P,M,AO,TV,TF,r,i,ni,R] = emboss(V,F,str,varargin)
  % EMBOSS Find a spot on a mesh to emboss with given text
  %
  % [VV,FF] = emboss(V,F,str)
  % [VV,FF,J] = emboss(V,F,str,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   str  string to emboss
  % Outputs:
  %   VV  #VV by 3 list of mesh vertex positions
  %   FF  #FF by 3 list of triangle indices into VV
  %   J  #FF list of indices into 1:size(F,1)+size(EF,1) where EF are the faces
  %     of the embossed letters: J<=size(F,1) indicates an original birth face
  %

  th = normrow(max(V)-min(V))/(100/3);
  fontname = [];
  engrave = false;
  max_r = inf;
  black_list = [];

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'BlackList','Engrave','Height','FontName','MaxRadius'}, ...
    {'black_list','engrave','th','fontname','max_r'});
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

  if islogical(black_list)
    black_list = find(black_list);
  end
  FW = F(setdiff(1:end,black_list),:);

  % Generate 10000 uniformly random samples 
  [P,I] = random_points_on_mesh(V,FW,10000,'Color','blue');
  % Compute the normals per-face
  NW = normalizerow(normals(V,FW));
  % Find the 160 closest points
  [K,D] = knnsearch(P,P,'K',160);
  % Compute normals at each point
  NWIK = reshape(NW(I(K),:),[size(K) 3]);
  % Compute a weighted fall-off per point
  W = exp(-(4*D./max(D,[],2)).^2);
  % Compute the ambient occlusion per point
  AO = ambient_occlusion(V,F,P,NW(I,:),1000);
  M = sum(W.*acos(min(max(sum(NWIK.*permute(NW(I,:),[1 3 2]),3),-1),1)).^2,2)+AO;
  [r,i] = max(D(:,end).*(M==min(M)));
  r = min(r,max_r);

  ni = NW(I(i),:);
  [w,a] = axisanglebetween(ni,[0 0 1],[1 0 0]);
  R = axisangle2matrix(w,a);
  [TV,TF] = text_to_mesh(str,'FontName',fontname,'TriangleFlags',' ');
  TV = TV-0.5*(max(TV)-min(TV));
  TV = TV/max(normrow(TV));
  TV(:,3) = TV(:,3)/max(TV(:,3));
  TV = TV*diag([r r th])*R+P(i,:);

  %clf;
  %hold on;
  %tsurf(F,V,fsoft,'FaceColor',orange,'EdgeColor','none','FaceAlpha',0.5);
  %tsurf(TF,TV,fsoft,'FaceColor',blue,'EdgeColor','none','FaceAlpha',0.5);
  %quiver3(P(i,1),P(i,2),P(i,3),N(I(i),1),N(I(i),2),N(I(i),3));
  %hold off;
  %camlight;
  %axis equal;

  if engrave
    [VV,FF,J] = mesh_boolean(V,F,TV,TF,'minus');
  else
    [VV,FF,J] = mesh_boolean(V,F,TV,TF,'union');
  end

end
