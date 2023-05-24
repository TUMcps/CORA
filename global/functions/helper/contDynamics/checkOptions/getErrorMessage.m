function msg = getErrorMessage(id)
% getErrorMessage - returns the error message corresponding to id
%
% Syntax:
%    msg = getErrorMessage(id)
%
% Inputs:
%    id - identifier for error message
%
% Outputs:
%    msg - error message
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author:       Mark Wetzlinger
% Written:      26-May-2021
% Last update:  19-May-2023 (MW, fix for unit test)
% Last revision:---

%------------- BEGIN CODE --------------

% input preprocessing
inputArgsCheck({{id,'att','char'}});

% enable access to codex
global codex;
if isempty(codex)
    initErrorCodex;
end

% check if id in codex
ididx = find(ismember(codex.id,id));

% empty message if id cannot be found
if isempty(ididx)
    msg = '';
else
    msg = codex.text{ididx};
end

end

%------------- END OF CODE --------------