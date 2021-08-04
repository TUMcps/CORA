function res = test_interval_isIntersecting
% test_interval_isIntersecting - unit test function of isIntersecting
%    note: only interval-to-interval tested!
%
% Syntax:  
%    res = test_interval_isIntersecting
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


% 1. Empty case: isIntersecting has to be false
res_empty = true;
I_empty = interval();
I_fullD = interval(-rand(3,1),rand(3,1));

if isIntersecting(I_fullD,I_empty)
    res_empty = false;
end

% dimension mismatch
res_mismatch = true;
I1 = interval(-1,1);
I2 = interval(-rand(2,1),rand(2,1));
try
    isIntersecting(I1,I2);
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res_mismatch = false;
    end
end

% combine results
res = res_empty && res_mismatch;

if res
    disp('test_isIntersecting successful');
else
    disp('test_isIntersecting failed');
end

%------------- END OF CODE --------------


