function res = test_interval_norm
% test_interval_norm - unit test function of norm
%
% Syntax:  
%    res = test_interval_norm
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
% Written:      27-July-2021
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% init interval, where upper bound vertex defines max-norm
lb = [-2; -1];
ub = [3; 4];
I = interval(lb,ub);

% compute norms: 1, 2, Inf
normI_1 = norm(I,1);
normI_2 = norm(I,2);
normI_Inf = norm(I,Inf);

% check with correct result
res = normI_1 == 7 && normI_2 == 5 && normI_Inf == 4;


if res
    disp('test_norm successful');
else
    disp('test_norm failed');
end

%------------- END OF CODE --------------
