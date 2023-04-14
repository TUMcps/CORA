function res = testLongDuration_conHyperplane_contains
% testLongDuration_conHyperplane_contains - unit test function of contains
%
% Syntax:  
%    res = testLongDuration_conHyperplane_contains
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
% See also: -

% Author:       Victor Gassmann
% Written:      19-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% combine tests
[~,res] = evalc('testLongDuration_mptPolytope_contains');

%------------- END OF CODE --------------