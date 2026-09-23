function pbpaste_into_workspace(dryrun)
  if nargin < 1
    dryrun = false;
  end

  [status, txt] = system('pbpaste');
  if status ~= 0
    error('Failed to read clipboard using pbpaste.');
  end

  txt = strtrim(txt);
  if isempty(txt)
    error('Clipboard is empty.');
  end

  fprintf('Clipboard length: %d chars\n', strlength(string(txt)));

  preview_len = 20;
  preview = txt;
  if strlength(string(preview)) > preview_len
    preview = extractBefore(string(preview), preview_len + 1) + " ...";
  end
  fprintf('Preview:\n%s\n', preview);

  if dryrun
    fprintf('[DRY RUN] Nothing evaluated.\n');
    return;
  end

  evalin('base', txt);
end
