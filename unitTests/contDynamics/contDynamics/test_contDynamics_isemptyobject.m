function res = test_contDynamics_isemptyobject
% test_contDynamics_isemptyobject - unit test for emptiness check
%
% Syntax:
%    res = test_contDynamics_isemptyobject
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
assert(isemptyobject(contDynamics()));

% only name not enough
assert(isemptyobject(contDynamics('sys')));

% dimension given, but zero
assert(isemptyobject(contDynamics('sys',0)));

% non-zero dimension given
assert(~isemptyobject(contDynamics('sys',1)));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
