function list = add2list(list,name,status,checkfuncs,errmsgs,varargin)
% add2list - adds an element consisting of a name, its status, and
%    its validation functions to a params or options list
%
% Syntax:
%    list = add2list(list,name,status,checkfuncs)
%    list = add2list(list,name,status,checkfuncs,condfunc)
%
% Inputs:
%    list - list of parameters or options
%    name - name of parameter or option
%    status - 'mandatory', 'optional', 'default'
%    checkfuncs - function handle or cell array of function handles,
%                  for validation of parameter or option
%    errmsgs - cell array of char array (error messages for failed check)
%    condfunc - function handle in cell, conditional function which
%                decides whether given option is required / redundant
%
% Outputs:
%    list - updated list
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      25-January-2021
% Last update:  28-January-2021
%               02-July-2021 (enact errmsgs)
% Last revision:---

%------------- BEGIN CODE --------------

admissiblestatus = {'mandatory','optional','default'};

% check varargin
condfunc = {};
if ~isempty(varargin)
    condfunc = varargin{1};
end

% check types of name, status, and checkfuncs
if ~ischar(name)
    throw(CORAerror('CORA:configFile',...
        "Name for parameter or option has to be a char array."));
elseif length(strfind(name,'.')) > 1
    throw(CORAerror('CORA:configFile',...
        "Name for parameter or options can contain at most one '.'"));
elseif ~ismember(admissiblestatus,status)
    throw(CORAerror('CORA:configFile',...
        sprintf(['Status of parameter or options has to be: \n  ',...
        strjoin(admissiblestatus,',\n  '),'.'])));
elseif ~iscell(checkfuncs) || (~isempty(condfunc) && ~iscell(condfunc))
    throw(CORAerror('CORA:configFile',...
        "Validation functions have to be within a cell-array."));
elseif iscell(condfunc) && length(condfunc) > 1
    throw(CORAerror('CORA:configFile',...
        "There can only be one conditional validation function."));
elseif ~iscell(errmsgs)
    throw(CORAerror('CORA:configFile',...
        "Error messages have to be within a cell-array."));
end

% check correctness of error message identifiers
for i=1:length(errmsgs)
    if ~ischar(errmsgs{i})
        throw(CORAerror('CORA:configFile',...
            "Error Messages have to be char arrays."));
    elseif isempty(getErrorMessage(errmsgs{i}))
        % check if errmsgs{i} actually exists
        throw(CORAerror('CORA:configFile',...
            "Error Message Identifier does not exist."));
    end
end

% check correctness of input arguments for function handles
for i=1:length(checkfuncs)
    if ~isa(checkfuncs{i},'function_handle')
        throw(CORAerror('CORA:configFile',...
            "Validation functions have to be function handles."));
    elseif nargin(checkfuncs{i}) ~= 1
        throw(CORAerror('CORA:configFile',...
            "Validation functions can only have one input argument."));
    end
end
if ~isempty(condfunc)
    if ~isa(condfunc{1},'function_handle')
        throw(CORAerror('CORA:configFile',...
            "The conditional validation function has to be a function handle."));
    elseif nargin(condfunc{1}) ~= 0
        throw(CORAerror('CORA:configFile',...
            "The conditional validation function has to have no input arguments."));
    end
end


% addition to list
list.name{end+1,1} = name;
list.status{end+1,1} = status;
list.checkfuncs{end+1,1} = checkfuncs;
list.errmsgs{end+1,1} = errmsgs;
if ~isempty(condfunc)
    list.condfunc{end+1,1} = condfunc{1};
else
    list.condfunc{end+1,1} = {};
end


end

%------------- END OF CODE --------------
