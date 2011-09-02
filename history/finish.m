setpref('StartupDirectory','LastWorkingDirectory',pwd) 
try 
    preserve_history; 
catch EM 
    h=msgbox(EM.message,sprintf('Error: %s',EM.identifier),'error'); 
    uiwait(h); 
    return 
end
