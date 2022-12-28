function data = imdata(src)
% IMDATA  converts the image file content into data URI scheme (RFC-2397)
%
%   S = imdata('foo.png');      reads the image file content as bytes,
%                               and returns them as base64 encoded
%                               data scheme.
%
%   Result: 'data:image/png;base64,iVBORw0KGgoAAAAN...'
%
% Adapted from:
% Written 2018-11-05 by Thomas Spriestersbach (fdx1601@icloud.com)
%
    fmt = [];
    src=lower(src);
    if endsWith(src,'.png' )
        fmt = 'image/png';
    elseif endsWith(src,'.jpg' ) || endsWith(src,'.jpeg' )
        fmt = 'image/jpeg';
    elseif endsWith(src,'.bmp' )
        fmt = 'image/bmp';
    else
        error( 'Unsupported image format' );
    end

    fp = fopen(src,'r');
    B = fread(fp);
    fclose(fp);
    Bstr = base64encode(B);
    data = ['data:' fmt ';base64,' Bstr];
    %fid=fopen(src,'r');
    %bytes=fread(fid,'uint8=>uint8');
    %fclose(fid);
    %data_old = ['data:' fmt ';base64,' matlab.net.base64encode(bytes)];
    %assert(isequal(data_old, data))
end

