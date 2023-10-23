function checks = add2checks(checkfun, errorId)
% add2checks - constructs an additional entry in the check struct of
%    dynamic parameter
%
% Syntax:
%    checks = add2checks(checkfun, errorId)
%
% Inputs:
%    checkfun - check function handle
%    errorId - error identified
%
% Outputs:
%    checks - struct
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: checkDynParameter, getErrorMessage

% Authors:       Tobias Ladner
% Written:       09-October-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% check input

% checkfun
if ~isa(checkfun,'function_handle')
    throw(CORAerror('CORA:configFile', "Validation functions have to be function handles."));
elseif nargin(checkfun) ~= 1
    throw(CORAerror('CORA:configFile', "Validation functions can only have one input argument."));
end

% errorId
if ~ischar(errorId)
    throw(CORAerror('CORA:configFile', "Error Messages have to be char arrays."));
elseif isempty(getErrorMessage(errorId))
    % check if errorId} actually exists
    throw(CORAerror('CORA:configFile',"Error Message Identifier does not exist."));
end

% set struct
checks = struct;
checks.checkfun = checkfun;
checks.errorId = errorId;
    
end

% ------------------------------ END OF CODE ------------------------------
