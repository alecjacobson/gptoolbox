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
clear lwd;
format short g;
