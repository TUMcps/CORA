function res = test_interval_convHull
% test_interval_convHull - unit test function of convHull
%
% Syntax:
%    res = test_interval_convHull
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

% Author:       Mark Wetzlinger
% Written:      12-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


% 1. Empty case: convHull is empty set
res_empty = true;
I_empty = interval();
I_fullD = interval(-rand(3,1),rand(3,1));

if ~isequal(convHull(I_fullD,I_empty),I_empty)
    res_empty = false;
end

% dimension mismatch
res_mismatch = true;
I_dim3 = interval(-rand(3,1),rand(3,1));
I_dim7 = interval(-rand(7,1),rand(7,1));
try
    convHull(I_dim3,I_dim7);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res_mismatch = false;
    end
end

% combine results
res = res_empty && res_mismatch;

if res
    disp('test_convHull successful');
else
    disp('test_convHull failed');
end

%------------- END OF CODE --------------

