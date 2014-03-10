function line = eat_comments(file_pointer,comment_symbol)
  % EAT_COMMENTS use fscanf to eat lines starting with a comment_symbol
  % 
  % line = eat_comments(file_pointer,comment_symbol)
  %
  % Inputs:
  %   file_pointer  returned from fopen
  %   comment_symbol  symbol of comment, e.g. '#'
  % Output:
  %   line  next line that does not start with a comment_symbol 
  %     (file_pointer is correspondingly adjusted)
  %
  assert(size(comment_symbol,2) == 1);

  while(true)
    % read next whole line
    line = fscanf(file_pointer,' %[^\n]s');
    if(size(line,2) == 0)
      break;
    elseif(line(1,1) ~= comment_symbol)
      break;
    end
  end
end
