function V = base64decode(RES)
  % V = base64decode(RES)
  % 
  % Inputs:
  %   RES  #RES long list of characters
  % Outputs:
  %   V  #V long list of bytes
  %
  % Faster implementation of:
  %  Documentation for matlab.net.base64decode
  %
  % See also: base64encode
  vals = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=';
  equ = find(RES=='=')-1;
  if isempty(equ)
    equ = numel(RES);
  end
  [~,inds] = ismember(RES(1:equ),vals);
  bin = dec2bin_logical(inds-1,6);
  %assert(isequal(bin,dec2bin(inds-1,6)=='1'));

  bin = reshape(bin',1,[]);
  bin = bin(1:floor(numel(bin)/8)*8);
  bin = reshape(bin,8,[])';
  % This is very slow...
  %V = bin2dec(bin);
  V = (bin)*(2.^(7:-1:0))';
end

