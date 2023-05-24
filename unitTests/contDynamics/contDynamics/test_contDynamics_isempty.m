function res = test_contDynamics_isempty
% test_contDynamics_isempty - unit test for emptiness check
%
% Syntax:
%    res = test_contDynamics_isempty
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

% Author:       Mark Wetzlinger
% Written:      16-May-2023
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% empty object
res = isempty(contDynamics());

% only name not enough
res(end+1,1) = isempty(contDynamics('sys'));

% dimension given, but zero
res(end+1,1) = isempty(contDynamics('sys',0));

% non-zero dimension given
res(end+1,1) = ~isempty(contDynamics('sys',1));

% combine results
res = all(res);

%------------- END OF CODE --------------