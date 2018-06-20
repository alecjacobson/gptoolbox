function [tf,te,to,sc,I,up] = tsurf_cad(F,V,varargin)
  % [tf,te,to,sc,I,up] = tsurf_cad(F,V,varargin)
  %
  % TSURF_CAD Display a triangle mesh in a CAD-like rendering style: edges,
  % shadow, floor.
  %
  % Inputs:
  %   F  #F by 3 list of face indices
  %   V  #V by 3 list of vertex positions
  % Outputs:
  %   tf  plot handle to mesh 
  %   te  plot handle to sharp edges 
  %   to  plot handle to boundary edges 
  %   sc  plot handle to shadow
  %   I  indices into V to cut mesh along sharp edges
  %   up  callback function to reset lines based on rotation
  %

  function set_hold(v)
    if v
      hold on;
    else 
      hold off;
    end
  end
  old_hold = ishold;
  %V = V*axisangle2matrix([1 0 0],-pi/2);
  N = normals(V,F);
  BC = barycenter(V,F);

  E = sharp_edges(V,F);
  E = union(sort(E,2),sort(outline(F),2),'rows');

  % cut mesh at sharp edges to get crisp normals
  [G,I] = cut_edges(F,E);
  W = V(I,:);

  if ~old_hold
    clf;
  end
  hold on;
  tf = tsurf(G,W, ...
    'FaceVertexCData',repmat(blue,size(W,1),1), ...
    'SpecularStrength',0, ...
    'DiffuseStrength',0.1, ...
    'AmbientStrength',1.0, ...
    'EdgeColor','none','FaceAlpha',0.9,fphong, ...
    varargin{:});
  te = tsurf(E(:,[1 2 2]),V,'EdgeColor',blue*0.75);
  to = tsurf([1 1 1],V,'LineWidth',2,'EdgeColor',blue*0.5);
  view(-9,34);
  axis equal;
  l = light('Position',[3 -4 5.0],'Style','infinite');
  [h,~,M,g] = add_shadow(tf,l,'Ground',[0 0 -1 min(V(:,3))-2e-3],'Fade','local','Color',[0.8 0.8 0.8]);
  % faint amient occlusion
  AO = ambient_occlusion(W,G,W,per_vertex_normals(W,G),1000);
  AO = AO*0.27;
  tf.FaceVertexCData = bsxfun(@times,tf.FaceVertexCData,1-AO);
  %set_hold(old_hold);

  
  % Hack so that axis doesn't change if V contains unreferenced points and
  % doesn't change.
  [BB,BF] = bounding_box(V);
  hold on;
  bf = trisurf(BF,BB(:,1),BB(:,2),BB(:,3),'CData',nan*BB(:,1),'EdgeColor','none');
  set_hold(old_hold);

  % floor board
  SV = [V ones(size(V,1),1)]*M';
  SV = bsxfun(@rdivide,SV(:,1:3),SV(:,4));

  [BB,BF] = bounding_box(SV(:,1:2));
  M = mean(bounding_box(V(:,1:2)));
  BBmM = bsxfun(@minus,BB,M);
  BB = bsxfun(@plus,sign(BBmM)*max(abs(BBmM(:))),M);
  BB = bsxfun(@plus,bsxfun(@minus,BB,mean(BB))*1.1,mean(BB));
  BB(:,3) = min(V(:,3))-4e-3;
  BB = reshape(permute(BB,[1 3 2]),[2 2 3]);
  % checkboard texture
  nx = 16;
  extent = @(BB) max(max(BB)) - min(min(BB));
  ny = ceil(extent(BB(:,:,2))/ extent(BB(:,:,1))*nx/2)*2;
  ny(isinf(ny)) = nx;
  ch = repmat(0.99-0.15*xor((mod(repmat(0:8*nx-1,8*ny,1),8*2)>7), ...
    (mod(repmat((0:8*ny-1)',1,8*nx),8*2)>7)),[1 1 3])*0.5 + 0.5;
  hold on;
  sc = surf(BB(:,:,1),BB(:,:,2),BB(:,:,3), ...
    'CData',ch,'FaceColor','texturemap', ...
    'SpecularStrength',0, 'DiffuseStrength',0, 'AmbientStrength',1);
  set_hold(old_hold);

  axis vis3d;
  camproj('persp');
  set(gca,'Visible','off');
  T = get(gca,'tightinset');
  set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]); 

  % Set up rotatation callbacks to hide view-dependent effects during drag
  up = @() ...
    set(to,'Faces', ...
      outline(F((sum(N.*bsxfun(@minus,BC,campos),2)<=0),:))*[1 0 0;0 1 1]) | ...
    cell2mat(cellfun(@(h) set(h,'FaceAlpha',0.5*(g*[campos 1]'<0)),h,'UniformOutput',false)) | ...
    set(sc,'FaceAlpha',1.0*(g*[campos 1]'<0));
  up();
  down = @() ...
    (exist('to','var') && ishandle(to) && isempty(set(to,'Faces',[]))) || ...
    isempty(set(rotate3d,'ActionPostCallback',[],'ActionPreCallback',[]));
  set(rotate3d,'ActionPostCallback',@(src,obj) up());
  set(rotate3d,'ActionPreCallback',@(src,obj) down());

%   for t = linspace(0,-360,60)
%    view(64+t,20);
%    up();
%    drawnow;
%    filename = 'tsurf-cad.gif';
%    figgif(filename,'nodither');
%   end
end
