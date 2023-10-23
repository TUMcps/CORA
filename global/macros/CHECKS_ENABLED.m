function res = CHECKS_ENABLED()
% CHECKS_ENABLED - macro to enable/disable input argument checks (these
%    checks are useful to ensure correct functionality and error tracing,
%    but are time-consuming and therefore not required if the calling code
%    does not contain any errors)
%
% Syntax:
%    res = CHECKS_ENABLED()
%
% Inputs:
%    ---
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: inputArgsCheck.m, equalDimCheck.m

% Authors:       Mark Wetzlinger
% Written:       02-November-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% ------------------------------ END OF CODE ------------------------------
