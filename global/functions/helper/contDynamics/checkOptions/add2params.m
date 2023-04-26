function add2params(name,status,checkfuncs,errmsgs,varargin)
% add2params - wrapper for add2list function, adds element to temporary
%    global variable fullParamsList
%
% Syntax:
%    add2params(name,status,checkfuncs,errmsgs,varargin)
%
% Inputs:
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
% Last revision:---

%------------- BEGIN CODE --------------

global fullParamsList;

if isempty(varargin)
    fullParamsList = add2list(fullParamsList,name,status,checkfuncs,errmsgs);
else
    fullParamsList = add2list(fullParamsList,name,status,checkfuncs,errmsgs,varargin{1});
end

end

%------------- END OF CODE --------------
