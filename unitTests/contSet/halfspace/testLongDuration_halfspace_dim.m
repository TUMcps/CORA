function res = testLongDuration_halfspace_dim
% testLongDuration_halfspace_dim - unit test function of dim
%
% Syntax:  
%    res = testLongDuration_halfspace_dim
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
% Written:      27-Sep-2019
% Last update:  16-March-2021 (MW, add empty case)
% Last revision:---

%------------- BEGIN CODE --------------
% Random cases
nrTests = 1000;
res_rand = true;

for i=1:nrTests
    % random dimension
    n = randi(50);
    
    % random normal vector and distance
    c = randn(n,1);
    d = randn(1);
    
    % init halfspace
    h = halfspace(c,d);
    
    % check result
    if dim(h) ~= n
        res_rand = false; break;
    end
end


% combine tests
res = res_rand;

if res_rand
    disp('test_dim successful');
else
    disp('test_dim failed');
end

%------------- END OF CODE --------------