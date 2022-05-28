function [I,A,FI,B] = shader(orig_tsh)
  % Act like a fragment shader.
  %
  % [I,A,FI,B] = shader(orig_tsh)
  %
  % Input:
  %   orig_tsh  handle to trisurf (must be triangle mesh)
  % Outputs:
  %   I  height by width by 3 rgb image frame from gcf window
  %   A  height by width boolean image of where mesh is
  %   FI height by width image of indices into rows of orig_tsh.Faces
  %   B  height by width by 3 image of corresponding barycentric coordinates
  %   

  I = im2double(getfield(getframe(gcf),'cdata'));

  orig_gcf = gcf;
  tsh = copy(orig_tsh);
  tsh = copyobj(tsh,gca);
  orig_visibile = orig_tsh.Visible;
  orig_tsh.Visible = 'off';

  orig_graphics_smoothing = get(gcf,'GraphicsSmoothing');
  set(gcf,'GraphicsSmoothing','off');

  set(tsh, ...
    'FaceColor','flat', ...
     'FaceVertexCData',repmat([1 0 1], size(tsh.Faces,1),1));
  C = getfield(getframe(gcf),'cdata');
  A = ((C(:,:,1)==255) & (C(:,:,2)==0) & (C(:,:,3)==255) );

  % encode integer as color
  rgb2id = @(I) sum(uint32(uint8(I)).*(uint32(256).^uint32(cat(3,2,1,0))),3);
  %[i,j,k] = ind2sub([256 256 256],65537+1);cat(3,i,j,k)-1
  % oof, do it as a lookup because anonymous functions need to be one line and
  % can't use multiple output parameters of subroutines...
  [GI,GJ,GK] = meshgrid(uint8(0:255),uint8(0:255),uint8(0:255));
  id2rgb = @(I) cat(3,GK(I+1),GI(I+1),GJ(I+1));
  %assert(all(0:256^3-1 == rgb2id(id2rgb(0:256^3-1))))

  if nargout>=3
    % Grad id 
    set(tsh, ...
      'FaceColor','flat', ...
      'FaceVertexCData',squeeze(id2rgb((1:size(tsh.Faces,1))')));
    FI = getfield(getframe(gcf),'cdata');
    FI = rgb2id(FI);
    FI = FI.*A;

    if nargout>=4
      FF = tsh.Faces;
      VV = tsh.Vertices;
      set(tsh, ...
        'Faces',reshape(1:numel(FF),[],3),'Vertices',VV(FF,:), ...
        'FaceColor','interp','FaceVertexCData',repdiag(ones(size(FF,1),1),3));
      B = im2double(getfield(getframe(gcf),'cdata'));
      B = B ./ sum(B,3);
      B = B.*A;
    end
  end

  %f = figure('WindowSTyle','docked');
  %imshow(A);
  %pause
  %close(f);

  orig_tsh.Visible = orig_visibile;
  % Why can't I just do delete(tsh) ?  Mauybe I can...
  children = get(gca,'Children');delete(children(1));
  figure(orig_gcf);
  set(gcf,'GraphicsSmoothing',orig_graphics_smoothing);
end
