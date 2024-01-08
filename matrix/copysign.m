function C = copysign(A,B)
  % C = copysign(A,B)
  %
  % Inputs:
  %   A  m by n matrix of doubles
  %   B  m by n matrix of doubles
  % Outputs:
  %   C  m by n matrix of doubles so that C(i) =  abs(A(i)) if B(i) >= 0 and
  %                                       C(i) = -abs(A(i)) if B(i)  < 0

  %C = (B>=0).*abs(A) + (B<0).*-abs(A);

  signmask = bitxor(typecast(double(0),'uint64'),typecast(double(-0),'uint64'));
  not_signmask = bitcmp(signmask,'uint64');

  %S = bitand(typecast(B,'uint64'),signmask);
  %A_abs = bitand(typecast(A,'uint64'), not_signmask);
  %C = typecast(bitor(A_abs,S),'double');

  % Combine into a single line
  C = typecast(bitor(bitand(typecast(A,'uint64'), not_signmask),bitand(typecast(B,'uint64'),signmask)),'double');
end
