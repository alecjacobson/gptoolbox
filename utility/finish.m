% Custom finish.m which is called on `quit`, saves history and current
% workspace.
%
% See also: startup.m

% Use preferences
% [instead](http://www.mathworks.ch/ch/help/matlab/matlab_env/exit-matlab-.html#bs6j5on-4)
% Preferences > MATLAB > General > Confirmation Dialogs > Confirm before
% exiting MATLAB
%
% if usejava('desktop')
%   button = questdlg('Are you sure you want to quit?', 'Exit Dialog','Cancel','Quit','Cancel');
%   
%   if strcmp(button,'Cancel') || strcmp(button,'')
%     quit cancel;
%     return
%   end
% end

lastworkspace = '/var/tmp/lastworkspace.mat';
if is_writable(lastworkspace)
  disp(['Saving workspace data to ' lastworkspace]);
  save(lastworkspace);
else
  warning('Workspace recovery location not writable');
end

if usejava('desktop')
  setpref('my','LastWorkingDirectory',pwd)
  
  if ~exist('preserve_history_failed_','var')
      preserve_history_failed_ = false;
  end
  try
      if ~preserve_history_failed_
          preserve_history;
      end
  catch
      preserve_history_failed_ = true;
      disp Exiting Matlab was halted so you can investigate.
      disp Type "preserve_history" after you fix it before exiting this Matlab session.
  end
end
