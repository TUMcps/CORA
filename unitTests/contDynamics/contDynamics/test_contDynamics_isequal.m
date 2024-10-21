function res = test_contDynamics_isequal
% test_contDynamics_isequal - unit test for equality check
%
% Syntax:
%    res = test_contDynamics_isequal
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
% Written:       19-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% empty object
assert(isequal(contDynamics(),contDynamics()));

% init contDynamics objects
sys1 = contDynamics('sys');
sys2 = contDynamics('sys',0);
sys3 = contDynamics('sys',1,1,1);

% compare systems
assert(isequal(sys1,sys2));
assert(~isequal(sys1,sys3));
assert(isequal(sys3,sys3));

% combine results
res = true;

% ------------------------------ END OF CODE ------------------------------
