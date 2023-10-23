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
res = isemptyobject(contDynamics());

% only name not enough
res(end+1,1) = isemptyobject(contDynamics('sys'));

% dimension given, but zero
res(end+1,1) = isemptyobject(contDynamics('sys',0));

% non-zero dimension given
res(end+1,1) = ~isemptyobject(contDynamics('sys',1));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
