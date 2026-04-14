function vobj = figmp4(name,vobj,nf)
  % vobj = figmp4(name,vobj)
  %
  % Inupts:
  %   name  path to .mp4 file
  %   vobj  [] on first call and see output for subsequent calls
  % Outputs:
  %   vobj  mp4 video object. the video will not be valid until you call
  %     vobj.close
  % 
  % Edit your code in ** three ** places
  %
  %   vobj = [];
  %   for iter = 1:10
  %      ...
  %      drawnow;
  %      vobj = figmp4('myvideo.mp4',vobj);
  %   end
  %   vobj.close;
  %
  if ~exist('nf','var')
    nf = 1;
  end
  if isempty(vobj)
    vobj = VideoWriter(name,'MPEG-4');
    vobj.Quality = 100;
    %vobj = VideoWriter(name,'Archival');
    vobj.open;
  end
  im = getfield(getframe(gcf),'cdata');
  % crop so that size(im,1) and size(im,2) are a multiple of 4
  sm = 16;
  sz = sm*floor([size(im,1) size(im,2)]/sm);
  im = im(1:sz(1),1:sz(2),:);

  for f = 1:nf
    vobj.writeVideo(im);
  end
end
