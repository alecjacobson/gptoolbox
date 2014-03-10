function definitions = tabcomplete(funcName, varargin)
%tabcomplete Sets tab completion for the specified function
%
% Syntax:
%    definitions = tabcomplete(funcName, argType, ...)
%
% Description:
%    Argument tab-completion occur when the user presses the <Tab> key
%    following the function name in the Command Window. The list of
%    available/possible arguments is then presented in a popup window.
%
%    TABCOMPLETE modifies the [matlabroot '/toolbox/local/TC.xml'] file
%    to set/unset the tab-completion definition of the specified function.
%    Note: changes take effect only after a Matlab restart.
%
%    TABCOMPLETE(funcName, argType1, argType2, ...) sets the tab-completion
%    list for function arguments 1, 2, etc. argType can be one of these pre-
%    defined keywords:
%      - 'var'       - list of current workspace variables
%      - 'fun'       - list of accessible functions
%      - 'subfun'    - list of accessible sub-functions
%      - 'dir'       - list of accessible folders (directories)
%      - 'file'      - list of accessible files (of any type)
%      - 'mfile'     - list of accessible *.m files (Matlab functions)
%      - 'matfile'   - list of accessible *.mat files (Matlab data)
%      - 'figfile'   - list of accessible *.fig files (figures)
%      - 'mdlfile'   - list of accessible *.mdl files (Simulink models)
%      - 'mcospkg'   - list of accessible MCOS class packages (R2010a+)
%      - 'mcosclass' - list of accessible MCOS classes (R2010a+)
%      - 'messageid' - list of accessible error/warning message IDs (R2011b+)
%      - (if none of the above is specified, 'subfun' is set automatically)
%    And in addition to the pre-defined kewords:
%      - additional string value(s) that will be added to the pop-up list
%
%    Multiple types for the same argument can be specified:
%      - as a cell array             - example: {'var','file','on','off'}
%      - as a space-separated string - example: 'var file on off'
%
%    The last argType specified will apply for all subsequent arguments.
%    The last argType may be '' (empty string) to indicate the end of tab-
%    completed args. The last argType does *NOT* accept additional string
%    values like the rest of the arguments. 
%
%    TABCOMPLETE(funcName,'') removes all arg tab-completions for funcName
%    (this is an immediate corollary of the previous paragraph).
%
%    definitions=TABCOMPLETE (without any arguments) returns a structure
%    array listing all the current tab-completion definitions.
%
%    definitions=TABCOMPLETE(funcName) returns a structure listing the
%    tab-completion definition for the specified function (if existing).
%
%    definitions=TABCOMPLETE(...) returns a structure array listing all
%    tab-completion definitions BEFORE their modification by TABCOMPLETE.
%
% Examples:
%    tabcomplete('addpath','dir')   % sets multiple folder completions
%    tabcomplete('edit','file fun') % sets multiple file/func completions
%    tabcomplete('cd','dir','')     % sets a single-arg folder completion
%    tabcomplete('save','matfile',{'var','-struct','-regexp'})
%       => *.mat followed by multiple variable-names or '-struct'/'-regexp'
%    tabcomplete('myFunc','')       % removes myFunc's arg tab-completions
%    defs = tabcomplete('cd')       % returns 'cd' function's tab-completions
%    defs = tabcomplete       % returns all currently-defined tab-completions
%
% Known issues/limitations:
%    - The modified tab-completions only take effect after Matlab restart.
%    - The last (default) argType does NOT accept non-pre-defined keywords
%      (this is a limitation of Matlab's TC.xsd definition file).
%    - Arguments *MUST* use at least one of the pre-define keywords for their
%      argType (this is another limitation of Matlab's TC.xsd definition file).
%    - Only lowercase function names are supported. This is another limitation
%      of Matlab's TC.xsd definition file. It can be overcome by editing TC.xsd
%      (in the same folder as TC.xml) line #20: Change <xsd:pattern value=
%      '[a-z_0-9]+(/[a-z_0-9]+)?'/> to: '[A-Za-z_0-9]+(/[A-Za-z_0-9]+)?'/>
%    - This utility currently does not enable setting previous-arg-based
%      tab-completions (e.g., as for 'whos' which defines: previous="-file").
%      You can edit the TC.xml file manually to achieve this functionality.
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
% Technical explanation:
%    A technical explanation of the code in this utility can be found on
%    <a href="http://undocumentedmatlab.com/blog/setting-desktop-tab-completions/">http://undocumentedmatlab.com/blog/setting-desktop-tab-completions</a>
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Change log:
%    2010-03-03: First version posted on the <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%    2010-09-28: Minor fix to msgbox in case of backup problem
%    2012-02-04: Added MCOSPKG,MCOSCLASS (R2010a+), MESSAGEID (R2011b+) argTypes; warn when using extra vars in R2010a-R2011b; attempt to make TC.xml writeable before updating
%
% See also:
%    cprintf, setPrompt (on the File Exchange)

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.02 $  $Date: 2012/02/04 19:33:47 $

    % Get the current definitions
    [defs,tcXmlFilename] = getCurrentXmlDefs;

    % If modification was requested
    if nargin > 1
        % Ensure that a backup file was made
        bakFilename = strrep(tcXmlFilename,'.xml','.bak');
        if ~exist(bakFilename,'file') && ~ispref(mfilename,'dontBackup')
            % Ask the user whether to save a backup copy (YES, no, no & don't ask again)
            msg = {'It is advisable to save a backup copy before modifying.', '', ['Save ' bakFilename '?']};
            switch getQuestDlg(msg)
                case 'Yes'   % => Yes: backup file
                    [status,msg] = copyfile(tcXmlFilename,bakFilename,'f');
                    if ~status
                        error('YMA:TabComplete:CreateBakup','Failed to create %s:\n%s',bakFilename,msg);
                    end
                case 'No & never ask again'   % => No & don't ask again
                    setpref(mfilename,'dontBackup',1);
                case ''
                    return;  % bail-out upon <Esc>
                otherwise
                    % forget it...
            end
        end
        
        % Find the relevant function index in the TC list
        funcIdx = find(strcmpi({defs.functionName},funcName));
        if length(funcIdx) > 1
            error('YMA:TabComplete:TooManyPossibilities','Too many possible functions - bailing out');
        elseif ~isempty(funcIdx) && ~isempty(defs(funcIdx).functionArgs)
            funcIdx = -funcIdx;  % indicate that there are <arg> sub-nodes
        end

        % Parse the requested changes and convert into an XML node
        nodeXmlStr = getXmlStr(funcName, varargin{:});

        % Update the XML file with the modified node
        fullXmlStr = getFullXmlStr(tcXmlFilename, funcName, funcIdx, nodeXmlStr);
        updateXmlFile(tcXmlFilename,fullXmlStr);
    end

    % Return the definitions if requested
    %if nargout
        if nargin == 1  % single arg: funcName
            definitions = defs(strcmpi({defs.functionName},funcName));
        else  % 0 or 2+ args: return entire defs list
            definitions = defs;
        end
    %end  % if nargout
end  % tabcomplete

% Query dialog
function answer = getQuestDlg(msg)
    createStruct.Interpreter = 'none';
    createStruct.Default = 'Yes';
    answer = questdlg(msg,mfilename,'Yes','No','No & never ask again',createStruct);
    drawnow;
end  % getQuestDlg

% Get the current definitions
function [defs,tcXmlFilename] = getCurrentXmlDefs

    % Check whether the TC.xml file exists
    tcXmlFilename = fullfile(matlabroot,'/toolbox/local/TC.xml');
    if ~exist(tcXmlFilename,'file')
        error('YMA:TabComplete:FileNotFound',strrep([tcXmlFilename ' was not found - cannot set tab-completions'],'\','\\'));
    end

    % Try parsing the defs file as an XML file
    try
        xmlDoc = xmlread(tcXmlFilename);
    catch
        error('YMA:TabComplete:XML_Error',strrep([tcXmlFilename ' appears to be an invalid XML file - parse error: ' lasterr],'\','\\'));  %#ok compatibility
    end

    % Loop over all child nodes
    try
        % Get the list of all TC binding elements.
        allBindings = xmlDoc.getElementsByTagName('binding');

        % Note that the item list index is zero-based.
        for k = 0 : allBindings.getLength-1

            % Get the base node's attributes
            thisBinding = allBindings.item(k);
            defs(k+1).functionName = getAttrValue(thisBinding, 'name');     %#ok grow
            defs(k+1).defaultType  = getAttrValue(thisBinding, 'ctype');    %#ok grow
            defs(k+1).extraValues  = getAttrValue(thisBinding, 'value');    %#ok grow
            defs(k+1).platform     = getAttrValue(thisBinding, 'platform'); %#ok grow
            defs(k+1).functionArgs = [];  %#ok grow

            % Loop over all binding args
            attrIdx = 1;
            childNode = thisBinding.getFirstChild;
            while ~isempty(childNode)
                %Filter out text, comments, and processing instructions.
                if childNode.getNodeType==childNode.ELEMENT_NODE && strcmpi(childNode.getNodeName,'arg')
                    argsIdx = str2num(getAttrValue(childNode, 'argn'));  %#ok
                    if isempty(argsIdx)
                        argsIdx = length(defs(k+1).functionArgs) + 1;
                        if ~isempty(getAttrValue(childNode, 'previous'))
                            argsIdx = max(2,argsIdx);
                        end
                    end
                    [defs(k+1).functionArgs(argsIdx).previousArg] = deal(getAttrValue(childNode, 'previous')); %#ok grow
                    [defs(k+1).functionArgs(argsIdx).argType]     = deal(getAttrValue(childNode, 'ctype'));    %#ok grow
                    [defs(k+1).functionArgs(argsIdx).extraValues] = deal(getAttrValue(childNode, 'value'));    %#ok grow
                    attrIdx = attrIdx + 1;
                end
                childNode = childNode.getNextSibling;
            end  % loop over all ars
        end  % loop over all bindings
    catch
        error('YMA:TabComplete:XML_Error','%s appears to be an invalid XML file - parse error:\n%s',strrep(tcXmlFilename,'\','\\'),lasterr);  %#ok compatibility
    end
end  % getCurrentXmlDefs

% Parse an XML attribute to get its value
function value = getAttrValue(node,name)
    attrStr = char(node.getAttributeNode(name));
    value = regexprep(attrStr,'.*="(.*)"','$1');
end  % getAttrValue

% Parse requested changes and convert into an XML node
function nodeXmlStr = getXmlStr(funcName, varargin)
    warnFlag = 0;
    if ~isempty(varargin{end})
        [ctype,extra,warnFlag] = getCtypeExtra(varargin{end},warnFlag,nargin-1); %#ok<ASGLU> extra is unused
        nodeXmlStr = ['  <binding name="' funcName '" ctype="' ctype '">'];
    else
        nodeXmlStr = ['  <binding name="' funcName '">'];
    end
    for argIdx = 1 : nargin-2
        [ctype,extra,warnFlag] = getCtypeExtra(varargin{argIdx},warnFlag,argIdx);
        if argIdx == 1,  nodeXmlStr = sprintf('%s\n',nodeXmlStr);  end
        nodeXmlStr = sprintf('%s    <arg argn="%d" ctype="%s"',nodeXmlStr,argIdx,ctype);
        if ~isempty(extra)
            V = sscanf(version, '%d.', 2);
            if V(1) == 7 && (V(2) >= 10 && V(2) <= 13) && ... % R2010a-R2011b
               ~ispref(mfilename,'dontUseExtra')
                msg = {['Matlab R2010a-R2011b have a bug that cause an endless loop (100% CPU load) when using extra tab-completion values such as "' extra '"'], ...
                       '', 'Use the extra parameters anyway?'};
                switch getQuestDlg(msg)
                    case 'Yes'   % => Yes: use extra, nothing to do here
                    case 'No & never ask again'   % => No & don't ask again
                        setpref(mfilename,'dontUseExtra',1);
                    otherwise    % Esc or No
                        extra = '';
                end
            end
        end
        if ~isempty(extra)
            nodeXmlStr = sprintf('%s value="%s"',nodeXmlStr,extra);
        end
        nodeXmlStr = sprintf('%s/>\n',nodeXmlStr);
    end  % loop over all args
    nodeXmlStr = sprintf('%s  </binding>',nodeXmlStr);
end  % getXmlStr

% Parse a values array to find CTYPE & extra values
function [ctype,extra,warnFlag] = getCtypeExtra(data,warnFlag,argIdx)
    values = normalizeData(data);
    if ~isempty(values)
        values = textscan(values,'%s ');
        values = unique(values{1});
    else
        values = {};
    end

    % Determine the accepted argTypes
    argTypes = {'VAR','FUN','SUBFUN','DIR','FILE','MFILE','MATFILE','FIGFILE','MDLFILE'};
    V = sscanf(version, '%d.', 2);
    if V(1) >= 8 || (V(1) >= 7 && V(2) >= 10)  % R2010a+
        argTypes = [argTypes, 'MCOSPKG', 'MCOSCLASS'];
    end
    if V(1) >= 8 || (V(1) >= 7 && V(2) >= 13)  % R2011b+
        argTypes = [argTypes, 'MESSAGEID'];
    end
    
    ctypeIdx = ismember(upper(values),argTypes);
    if ~any(ctypeIdx)
        if ~warnFlag
            warning('YMA:TabComplete:NoStdCtype','No standard type was specified for arg #%d - setting it to ''SUBFUN''', argIdx);
            warnFlag = 1;
        end
        ctype = 'SUBFUN';    % No standard ctype was specified so use default = 'SUBFUN'
    else
        ctype = normalizeData(upper(values(ctypeIdx)));
    end
    extra = normalizeData(values(~ctypeIdx));
end  % getCtypeExtra

% Normalize a cell array into a space-separated string
function str = normalizeData(str)
    try
        str = sprintf('%s ',str{:});
    catch
        % Never mind - probably not a cell array but a string
    end
    str = strtrim(str);
end  % normalizeData

% Update the XML text string with the modified node
function fullXmlStr = getFullXmlStr(tcXmlFilename, funcName, funcIdx, nodeXmlStr)
    txt = urlread(['file:///' tcXmlFilename]);  % simple hack to quickly read the file...
    if isempty(funcIdx)
        fullXmlStr = regexprep(txt,'(</TC>)',sprintf('%s\n$1',nodeXmlStr));
    else
        if funcIdx > 0
            searchStr = sprintf(' *<binding +name="%s"[^/]*/>',funcName);
            if isempty(regexp(txt,searchStr,'once'))
                searchStr = sprintf(' *<binding +name="%s".*?</binding>',funcName);
            end
        else
            searchStr = sprintf(' *<binding +name="%s".*?</binding>',funcName);
        end
        fullXmlStr = regexprep(txt,searchStr,sprintf('%s',nodeXmlStr));
    end
end  % getFullXmlStr

% Update the XML file with the modified node
function updateXmlFile(tcXmlFilename,fullXmlStr)
    try
        % First try to make the file writeable - try to proceed even upon error
        try fileattrib(tcXmlFilename,'+w'); catch, end
        [fid,msg] = fopen(tcXmlFilename,'wb');
        if fid<0,  error(msg);  end
        len = fprintf(fid, '%s', fullXmlStr);
        fclose(fid);
        if len<0,  error('Not all requested text was written');  end
    catch
        error('YMA:TabComplete:CannotSave','Could not update %s:\n%s',tcXmlFilename,lasterr);  %#ok compatibility
    end
end  % updateXmlFile

% Copyright (c) 2010, Yair Altman
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
