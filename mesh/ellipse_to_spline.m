function [Pi,Ci] = ellipse_to_spline(cx,cy,rx,ry)
  % special fraction
  s =  1193/2160;
  Pi = [1,0;1,s;s,1;0,1;-s,1;-1,s;-1,0;-1,-s;-s,-1;0,-1;s,-1;1,-s;1,0];
  % This repeats the first point rather than closing. parse_path may be counting
  % on this.
  Ci = [1,2,3,4;4,5,6,7;7,8,9,10;10,11,12,13];
  Pi = Pi.*[rx ry]+[cx cy];
end

