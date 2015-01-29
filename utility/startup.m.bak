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

if ispref('StartupDirectory','LastWorkingDirectory')
    lwd = getpref('StartupDirectory','LastWorkingDirectory');
    try
        cd(lwd)
    catch
        disp('Sorry, but I could not go to your last working directory:')
        disp(lwd)
    end;
end;
com.mathworks.mde.desk.MLDesktop.getInstance.restoreLayout('figure-command-history');
clear lwd;
format short g;
