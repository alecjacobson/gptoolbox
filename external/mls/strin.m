function contained = strin(str,strs,caseins)
% STRINI  Check if a string is contained in a set
%
%  This function checks if a string is contained in a set of strings (a
% cell). The returned value is a logical value.
%
%  Params
%  ------
% IN:
%  str          = The string to be found.
%  strs         = The cell of strings.
%  casein       = Use a case insensitive search? (def=false)
% OUT:
%  contained    = The string str is contained in strs?
%
%  Pre
%  ---
% -  The provided params must be all strings.
%
%  Post
%  ----
% -  The check i done and a correct result is returned without unhandled
%   exceptions (this to underline that if the pre aren't granted exceptions
%   can occour).
%
%  SeeAlso
%  -------
% strcmp, strcmpi
%
%  Examples
%  --------
%   Check a type of command:
% >> if strin(command,{'command1','command2'},true) ...

% Set default values
if nargin<3 caseins=false; end

% Select the comparing function
if caseins
    cmpfunc = @strcmpi;
else
    cmpfunc = @strcmp;
end

% Iterate on strings
for ind=1:numel(strs)
    % Check a single string
    if feval(cmpfunc,str,strs{ind})
        contained = true;
        return;
    end
end

% No one apply!
contained = false;
