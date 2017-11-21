% Custom startup which loads last working directory and workspace
%
% See also: finish

lastworkspace = '/var/tmp/lastworkspace.mat';
try
  load(lastworkspace);
catch
  disp('Sorry, but I could not load last workspace from:')
  disp(lastworkspace)
end;

if ispref('my','LastWorkingDirectory')
    lwd = getpref('my','LastWorkingDirectory');
    try
        cd(lwd)
    catch
        warning('Sorry, could not change to your last working directory: %s', lwd);
    end;
end;
clear lwd;

com.mathworks.mde.desk.MLDesktop.getInstance.restoreLayout('figure-command-history');
clear lwd;
format short g;
