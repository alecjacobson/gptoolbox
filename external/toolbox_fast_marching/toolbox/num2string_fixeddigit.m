function str = num2string_fixeddigit(num, d)

% num2string_fixeddigit - convert a number to string with a fixed number of digits
%
%   str = num2string_fixeddigit(num, d);
%
%   examples for d=3: 12.3->'012', 99->'099', 999->'999', 1200->'1200' (overflow).
%
%   Copyright (c) 2004 Gabriel Peyré


num = round(num);

str = num2str(num);

for i=1:d-1
    if num<10^i
        str = ['0' str];
    end
end