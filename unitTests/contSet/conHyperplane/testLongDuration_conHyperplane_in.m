function res = testLongDuration_conHyperplane_in
% testLongDuration_conHyperplane_in - unit test function of in
%
% Syntax:  
%    res = testLongDuration_conHyperplane_in
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
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
[~,res] = evalc('testLongDuration_mptPolytope_in');

if res
    disp('testLongDuration_conHyperplane_in successful');
else
    disp('testLongDuration_conHyperplane_in failed');
end

%------------- END OF CODE --------------