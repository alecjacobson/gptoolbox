button = questdlg('Are you sure you want to quit?', 'Exit Dialog','Yes','No','No');

if strcmp(button,'No')
  quit cancel;
  return
end

lastworkspace = '/var/tmp/lastworkspace.mat';
disp(['Saving workspace data to ' lastworkspace]);
save(lastworkspace);

setpref('StartupDirectory','LastWorkingDirectory',pwd) 
try 
    preserve_history; 
catch EM 
    h=msgbox(EM.message,sprintf('Error: %s',EM.identifier),'error'); 
    uiwait(h); 
    return 
end

