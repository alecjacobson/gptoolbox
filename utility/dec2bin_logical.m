function bin = dec2bin_logical(D,numBits)
  % bin = dec2bin_logical(D,numBits)
  %
  % Slightly faster version of dec2bin that outputs a logical array directly
  %
  % Documentation for dec2bin:
  %
  % dec2bin Convert decimal integer to its binary representation
  % dec2bin(D) returns the binary representation of D as a character
  % vector. D must be an integer. If D is greater than flintmax, dec2bin
  % might not return an exact representation of D.
  %
  % dec2bin(D,numBits) produces a binary representation with at least
  % numBits bits.
  %
  % Example
  %    dec2bin(23) returns '10111'
  % 
  % See also bin2dec, dec2hex, dec2base, flintmax.
  %
  % Example:
  %   D = [1;4;7];
  %   n = 3;
  %   bin = dec2bin_logical(D,n);
  %   dec = bin*2.^(n-1:-1:0)';
  %   assert(isequal(D,dec))

  D = reshape(D,[],1);
  assert(all(D<2^numBits),'numBits not large enough');

  bin = false(numel(D),numBits);
  for i = numBits:-1:1
    m = 2^(i-1);
    bin(:,numBits-i+1) = D>=m;
    D = D-bin(:,numBits-i+1)*m;
  end

end
