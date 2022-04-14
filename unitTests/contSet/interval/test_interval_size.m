function res = test_interval_size
% test_interval_size - unit test function of size
%
% Syntax:  
%    res = test_interval_size
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
% Written:      28-August-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% create interval
dim = 1 + floor(5 * rand(1));
temp1 = -rand(dim,1);
temp2 = rand(dim,1);
I = interval(min([temp1,temp2],[],2),max([temp1,temp2],[],2));

% convert to zonotope
[rows, cols] = size(I);

% true size
rows_true = dim;
cols_true = 1;

% compare results
res = rows == rows_true && cols == cols_true;


if res
    disp('test_size successful');
else
    disp('test_size failed');
end

%------------- END OF CODE --------------