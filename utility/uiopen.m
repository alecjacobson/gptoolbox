function uiopen(type,direct)
% UIOPEN overloaded for custom Files. Do not change the file name of this
% file. Remember you are overloading uiopen inside toolbox/matlab/uitools
%

%----Your file extension -----v
if     ((~isempty(findstr(type,'.off'))) && (direct))
elseif ((~isempty(findstr(type,'.obj'))) && (direct))
elseif ((~isempty(findstr(type,'.mesh'))) && (direct))
elseif ((~isempty(findstr(type,'.dmat'))) && (direct))
elseif ((~isempty(findstr(type,'.ply'))) && (direct))
elseif ((~isempty(findstr(type,'.stl'))) && (direct))
elseif ((~isempty(findstr(type,'.wrl'))) && (direct))
    %-------------------------------------------------
    % Your function that will open/run this file type
    %-------------------------------------------------

    %-------------------------------------------------
else
    %----------DO NOT CHANGE---------------------------
    presentPWD = pwd;
    cd([matlabroot '/toolbox/matlab/uitools']);
    strn = ['uiopen(''' type ''',' num2str(direct) ')'];
    eval(strn);
    cd(presentPWD);
    %----------DO NOT CHANGE---------------------------
end
%-------------------------------------------------------------------------
