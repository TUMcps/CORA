function res = test_interval_and
% test_interval_and - unit test function of logical conjunction,
%    overloaded '&' operator for intervals
%
% Syntax:  
%    res = test_interval_and
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
% See also: mtimes

% Author:       Dmitry Grebenyuk, Mark Wetzlinger
% Written:      05-January-2016
% Last update:  23-April-2023 (MW, add empty set cases)
% Last revision:---

%------------- BEGIN CODE --------------

res = true;

% 3D case
I1 = interval([-10; -2; 10], [-5; 8; 15]);
I2 = interval([-11; 2; 11], [-6; 9; 12]);
I = I1 & I2;
I_true = interval([-10;2;11],[-6;8;12]);
res(end+1,1) = isequal(I,I_true);

% empty intersection
I1 = interval([-5;-2],[2;4]);
I2 = interval([-7;6],[-3;8]);
I = I1 & I2;
res(end+1,1) = isempty(I);

% combine results 
res = all(res);

%------------- END OF CODE --------------