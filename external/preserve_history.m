function preserve_history ()
% Maintains command history indefinitely in a separate history_saved.m file
% since Matlab is truncating history.m to a fixed size. The code is fast even
% when history_saved.m gets very big.
%
% Also works with the new format introduced in R2014a which is XML. Its
% contents are parsed as plain commands and merged into history_saved.m
%
% Bogdan Roman, 2014, abr28@cam.ac.uk, University of Cambridge

savpath = fullfile(prefdir,'history_saved.m');

% If History.xml exists (R2014a or later) read that, else revert to history.m
hpath = fullfile(prefdir,'History.xml');
if exist(hpath,'file')
    % Read the XML file into a string to look as the older history.m (i.e.
    % remove all XML tags)
    mathist = fileread(hpath);
    mathist = regexprep(mathist, '(<[^>]+>\s*)+', '\n', 'lineanchors');
    % translate html entities and remove leading newline
    mathist = strrep(mathist(3:end), '&gt;', '>');
    mathist = strrep(mathist, '&lt;', '<');
else
    mathist = fileread(fullfile(prefdir,'history.m'));
end

% replace \r and \r\n with \n (safe to copy between OSes)
mathist = regexprep(mathist, '\r(\n)?', '\n');

fprintf('Merging current command history into %s ...\n', savpath);

% if empty saved file, dump it all
if ~exist(savpath,'file')
    fid = fopen(savpath, 'w');
    fwrite(fid, mathist);
    fclose(fid);
    return;
end

% Read saved file fully and use regexp(). Faster than parsing line by line.
savhist = fileread(savpath);
%savhist = regexprep(savhist, '\r(\n)?', '\n');

% Look for last saved date in the preserved file. Must be a valid date and
% not 'Unknown date', which may match in both files and overwrite the entire
% saved history. Note the \d to that effect.
[savstart, savend] = regexp(savhist,'^%-- \d.*?$','lineanchors');
if isempty(savstart)
    % If none found, then append current history entirely
    fid = fopen(savpath,'a');
    fwrite(fid, mathist);
    fclose(fid);
else
    % Extract last save date as string and search for it in current history
    savlastdate = savhist(savstart(end):savend(end)-1);
    matstart = regexp(mathist, ['^' savlastdate], 'lineanchors');
    % If not found, append the entire current history
    if isempty(matstart)
        fid = fopen(savpath,'a');
        fwrite(fid, mathist);
        fclose(fid);
    % Otherwise merge the chunks
    else
        fid = fopen(savpath, 'r+');
        fseek(fid, savstart(end)-1, 'bof');
        fwrite(fid, mathist(matstart(1):end));
        fclose(fid);
    end
end
disp('Done');
