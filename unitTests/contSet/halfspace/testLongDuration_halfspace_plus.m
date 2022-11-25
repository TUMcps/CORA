function res = testLongDuration_halfspace_plus
% testLongDuration_halfspace_plus - unit test function of plus
%
% Syntax:  
%    res = testLongDuration_halfspace_plus
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

% Random tests
res_rand = true;
nrTests = 1000;

for i=1:nrTests
    % random dimension
    n = randi(50);
    
    % random normal vector (unit length)
    c = randn(n,1);
    c = c / vecnorm(c,2);
    
    % random distance
    d = randn(1);
    
    % init halfspace
    h = halfspace(c,d);
    
    % random vector
    v = randn(n,1);
    
    % compute result
    h_plus = h + v;
    
    % true result
    d_true = d + c.' * v;
    h_true = halfspace(c,d_true);
    
    % compare results
    if ~isequal(h_plus,h_true)
        res_rand = false; break;
    end       
    
end


% combine tests
res = res_rand;

if res_rand
    disp('test_plus successful');
else
    disp('test_plus failed');
end

%------------- END OF CODE --------------