function S = last_command(varargin)
  % LAST_COMMAND returns the last issued command as a string.
  %
  % S = last_command();
  % S = last_command('ParamName',ParamValue,...);
  % 
  % Input (optional):
  %   'Offset' Followed by a number, rather than finding last command, finds
  %     last-offset command,{0}
  %   'LimitToSession' optionally followed by true or false, searches through
  %     all history {default} or just session history
  % Outputs:
  %   S  string of last command issued

  offset = 0;
  % look in all history (or just session)
  all_history = true;

  ii = 1;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Offset'
      ii = ii + 1;
      assert(ii<=nargin);
      offset = varargin{ii};
    case 'LimitToSession'
      if( (ii+1)<=nargin && ~ischar(varargin{ii+1}))
        ii = ii + 1;
        all_history = varargin{ii};
      else
        all_history = false;
      end
    otherwise
      error('Invalid parameter name');
    end
    ii = ii+1;
  end

  if all_history
    history = com.mathworks.mlservices.MLCommandHistoryServices.getAllHistory;
  else
    history = com.mathworks.mlservices.MLCommandHistoryServices.getSessionHistory;
  end
  S = history(end-offset);
end
