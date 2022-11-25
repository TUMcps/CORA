function res = test_halfspace_plus
% test_halfspace_plus - unit test function of plus
%
% Syntax:  
%    res = test_halfspace_plus
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
% Written:      16-March-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% 1. empty case
res_empty = true;
h = halfspace();
v = rand(2,1);
try
    h + v; % should throw error here
    res_empty = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:emptySet')
        res_empty = false;
    end
end


% 2. wrong input
res_input = true;
h = halfspace(randn(2,1),1);
v = rand(2,2);
try
    h + v; % should throw error here
    res_input = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:wrongValue')
        res_input = false;
    end
end


% 3. dimension mismatch
res_dim = true;
h = halfspace(randn(2,1),1);
v = rand(3,1);
try
    h + v; % should throw error here
    res_dim = false;
catch ME
    if ~strcmp(ME.identifier,'CORA:dimensionMismatch')
        res_dim = false;
    end
end



% combine tests
res = res_empty && res_input && res_dim;

if res
    disp('test_plus successful');
else
    disp('test_plus failed');
end

%------------- END OF CODE --------------