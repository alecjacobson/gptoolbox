function varargout = showvoronoi(D,V,E,RP,RD,crop)
  % SHOWVORONOI display a voronoi diagram
  %
  % h = showvoronoi(D,V,E,RP,RD)
  % h = showvoronoi(D,V,E,RP,RD,crop)
  % 
  % Inputs:
  %   D  #D by 2 list of voronoi cell centers (data points)
  %   V  #V by 2 list of voronoi diagram vertices
  %   E  #E by 2 list of voronoi internal edges
  %   RP  #RP by 1 list of voronoi infinite ray start points
  %   RD  #RP by 2 list of voronoi infinite direction vectors
  %   optional:
  %     crop  crop plot to fit just data points {true}
  % Output:
  %   h  handle to plot
  %
  % See also readEDGE, triangle
  %
  if(~exist('crop','var'))
    crop = true;
  end
  RD = normalizerow(RD);
  offset = 2*max(sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2)));
  EX = [V(E(:,1),1),V(E(:,2),1); V(RP,1) V(RP,1)+offset*RD(:,1)]';
  EY = [V(E(:,1),2),V(E(:,2),2); V(RP,2) V(RP,2)+offset*RD(:,2)]';
  h = plot(EX,EY,'-',D(:,1),D(:,2),'.');
  if(crop)
    set(h(1:end-1),'xliminclude','off','yliminclude','off');
  end
  if(nargout == 1)
    varargout{1} = h;
  end
end
