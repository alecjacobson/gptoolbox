function sqrD = segment_segment_squared_distance(S1,D1,S2,D2)
  % SEGMENT_SEGMENT_SQUARED_DISTANCE Compute the squared distance from the
  % points of closest approach on each pair of line segments: the Euclidean
  % distance between the line segments.
  %
  % Inputs:
  %   S1  #segments by dim list of source points on the "first" segments
  %   D1  #segments by dim list of destination points on the "first" segments
  %   S2  #segments by dim list of source points on the "second" segments
  %   D2  #segments by dim list of destination points on the "second" segments
  % Outputs:
  %   sqrD  #segments list of squared distances 
  %

  assert(size(S1,1) == size(D1,1));
  assert(size(S2,1) == size(D2,1));
  % Each list of endpoints is allowed to be a single segment: we'll just repeat
  % it to match the length of the other list.
  n1 = size(S1,1);
  n2 = size(S2,1);
  if n1 == 1
    S1 = repmat(S1,n2,1);
    D1 = repmat(D1,n2,1);
  end
  if n2 == 1
    S2 = repmat(S2,n1,1);
    D2 = repmat(D2,n1,1);
  end
  assert(size(S1,1) == size(S2,1));
  assert(size(D1,1) == size(D2,1));

  % Translated from:
  % http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment()
  u = D1 - S1;
  v = D2 - S2;
  w = S1 - S2;
  a = sum(u.*u,2);         % always >= 0
  b = sum(u.*v,2);
  c = sum(v.*v,2);         % always >= 0
  d = sum(u.*w,2);
  e = sum(v.*w,2);
  D = a.*c - b.*b;        % always >= 0
  sc = D;
  sN = D;
  sD = D;       % sc = sN / sD, default sD = D >= 0
  tc = D;
  tN = D;
  tD = D;       % tc = tN / tD, default tD = D >= 0

  % % compute the line parameters of the two closest points
  % if (D < SMALL_NUM) { % the lines are almost parallel
  %   sN = 0.0;         % force using point P0 on segment S1
  %   sD = 1.0;         % to prevent possible division by 0.0 later
  %   tN = e;
  %   tD = c;
  % } else {                 % get the closest points on the infinite lines
  %   sN = (b*e - c*d);
  %   tN = (a*e - b*d);
  %   if (sN < 0.0) {        % sc < 0 => the s=0 edge is visible
  %     sN = 0.0;
  %     tN = e;
  %     tD = c;
  %   } else if (sN > sD) {  % sc > 1  => the s=1 edge is visible
  %     sN = sD;
  %     tN = e + b;
  %     tD = c;
  %   }
  % }
  sN = b.*e - c.*d;
  tN = a.*e - b.*d;
  tN(sN<0) = e(sN<0);
  tD(sN<0) = c(sN<0);
  sN(sN<0) = 0;
  tN(sN>sD) =  e(sN>sD)+b(sN>sD);
  tD(sN>sD) =  c(sN>sD);
  sN(sN>sD) = sD(sN>sD);
  sN(D<=eps) = 0.0;
  sD(D<=eps) = 1.0;
  tN(D<=eps) = e(D<=eps);
  tD(D<=eps) = c(D<=eps);

  % if (tN < 0.0) {            % tc < 0 => the t=0 edge is visible
  %   tN = 0.0;
  %   % recompute sc for this edge
  %   if (-d < 0.0)
  %     sN = 0.0;
  %   else if (-d > a)
  %     sN = sD;
  %   else {
  %     sN = -d;
  %     sD = a;
  %   }
  % } else if (tN > tD) {      % tc > 1  => the t=1 edge is visible
  %   tN = tD;
  %   % recompute sc for this edge
  %   if ((-d + b) < 0.0)
  %     sN = 0;
  %   else if ((-d + b) > a)
  %     sN = sD;
  %   else {
  %     sN = (-d +  b);
  %     sD = a;
  %   }
  % }
  sN( tN<0 ) = -d(tN<0);
  sD( tN<0 ) =  a(tN<0);
  sN( (tN<0) & (-d<0) ) = 0;
  sN( (tN<0) & (-d>a) ) = sD( (tN<0) & (-d>a) );
  tN(tN<0) = 0;
  sN(tN>tD) = (-d(tN>tD) + b(tN>tD));
  sD(tN>tD) = a(tN>tD);
  sN( (tN>tD) & ((-d+b)<0) ) = 0;
  sN( (tN>tD) & ((-d+b)>a) ) = sD( (tN>tD) & ((-d+b)>a) );
  tN(tN > tD) = tD(tN > tD);

  % finally do the division to get sc and tc
  sc = sN./sD;
  sc(abs(sN)<=eps) = 0;
  tc = tN./tD;
  tc(abs(tN)<=eps) = 0;

  % get the difference of the two closest points
  % =  S1(sc) - S2(tc)
  dP = w + bsxfun(@times,sc, u) - bsxfun(@times,tc,v);

  sqrD = sum(dP.*dP,2);
end
