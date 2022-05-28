function RES = base64encode(V)
% RES = base64encode(V)
%
% Faster implementation of:
%  Documentation for matlab.net.base64encode
%
%  Perform Base 64 encoding of a string or vector of bytes
%    RES = base64encode(V) encodes a string, character vector, or numeric vector using
%    Base 64 encoding as documented in RFC 4648, section 4 and returns the encoded
%    characters as a string.  This encoding is used in a number of contexts in
%    Internet messages where data must be transmitted in a limited set of ASCII
%    characters.  It is often used to encode strings which may have special characters
%    that might be misinterpreted as control characters by the transmission protocol,
%    but it is also used to encode arbitrary binary data.
%
%    If the input is a string or character vector, it is first converted to bytes
%    using the user default encoding.  If you want to use a different character
%    encoding, use unicode2native to convert the string to a uint8 vector before
%    passing it into this function.
%
%  See also base64decode, unicode2native
%

  %% This is really slow and doing a stupid for loop:
  %Vstr = matlab.net.base64encode(V);
  vals = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=';
  bin = dec2bin_logical(V,8);
  %assert(isequal(bin, dec2bin(V,8)=='1'));

  bin = reshape(bin',1,[]);
  Vpad = mod(6-mod(numel(bin),6),6);
  bin = [bin false(1,Vpad)];
  inds = (2.^(5:-1:0))*(reshape(bin,6,[]))+1;
  Vstr = vals(inds);
  Vstr_pad = mod(4-mod(numel(Vstr),4),4);
  RES = [Vstr repmat('=',1,Vstr_pad)];
  %Vdec = base64decode(RES);
  %assert(isequal(V,Vdec))

end
