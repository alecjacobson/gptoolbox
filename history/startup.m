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