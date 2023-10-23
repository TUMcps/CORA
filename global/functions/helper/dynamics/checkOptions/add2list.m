function S = add2list(name,status)
% add2list - constructs an additional entry in a params or options list
%    consisting of a name and its status. Default values etc. are then 
%    set centrally, see 'See also' below.
%
% Syntax:
%    S = add2list(name,status)
%
% Inputs:
%    name - name of parameter or option
%    status - 'mandatory', 'optional', 'default'
%
% Outputs:
%    S - struct with name, status, checkfuncs, errmsgs, condfunc
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: getDefaultValue,checkDynParameter, getMembers, getErrorMessage,
%     getCondfunDynParameter, postProcessing

% Authors:       Mark Wetzlinger
% Written:       25-January-2021
% Last update:   28-January-2021
%                02-July-2021 (MW, enact errmsgs)
%                19-June-2023 (MW, skip error messages, return struct)
%                05-October-2023 (TL, removed checksfuncs, errmsgs, condfunc)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

admissiblestatus = {'mandatory','optional','default'};

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
end

% init struct
S = struct();
S.name = name;
S.status = status;
S.condfun = []; % will be set later

% ------------------------------ END OF CODE ------------------------------
