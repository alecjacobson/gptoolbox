% JS - 11/01/2005
% Maintains command history indefinitely in file "history_saved.m"
%
% This is handy because Matlab limits the history file to 20k in size, so
% throws away old commands (which are invariably the ones you want to recall).
% To automate history preservation call this script from "startup.m" or
% "finish.m"

% find last date matlab was opened in 'preserved' history file
[fid, err] = fopen(fullfile(prefdir,'history_saved.m'), 'r'); 

i=0;
while ~feof(fid)
    i=i+1;
    line = fgetl(fid);
    if regexp(line,'%-- ') == 1
        date_line = line;
        date_line_pos = ftell(fid);
    end
end

%date_line;
fclose(fid);

% match last date from 'preserved' file to a line in current 'history.m' file, then append new lines in 'history.m' to 'preserved' file
[fid, err] = fopen(fullfile(prefdir,'history.m'), 'r'); 

found = 0;
i=0;
while ~feof(fid)
    i=i+1;
    line = fgetl(fid);
    whole_file{i,:} = line;
    if regexp(line,date_line) == 1
        alpha = 1;
        while ~feof(fid)
            found = 1;
            line = fgetl(fid);
            append{alpha,:} = line;
            alpha = alpha + 1;
        end
    end
end

fclose(fid);

if found == 1                   % add commands in 'append' to end of saved
    [fid, err] = fopen(fullfile(prefdir,'history_saved.m'), 'r'); 

    fseek(fid,date_line_pos,'bof');             % move to last recognised date in saved file
    for i=1:length(append)
        fprintf(fid, '%s\r\n', append{i});
    end

    fclose(fid);

    % elseif found == 0               % add ALL commands from current to end of saved
    %     fseek(fid,0,'eof');             % move to end of file
    %     for i=1:length(whole_file)
    %         fprintf(fid, '%s\r\n', whole_file{i});
    %     end
    %
end


disp('History Preserved');
