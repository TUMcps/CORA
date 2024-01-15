function res = test_intervalMatrix_interval
% test_intervalMatrix_interval - unit test function for conversion to
%    interval objects
% 
% Syntax:
%    res = test_intervalMatrix_interval
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

% Authors:       Mark Wetzlinger
% Written:       18-June-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true(0);

% only one row
c = [1 2 3 4];
d = [1 1 1 1];
intMat = intervalMatrix(c,d);
I = interval(intMat);
I_ = interval(c-d,c+d);
res(end+1,1) = isequal(I,I_);

% only one column
c = [1 2 3 4]';
d = [1 1 1 1]';
intMat = intervalMatrix(c,d);
I = interval(intMat);
I_ = interval(c-d,c+d);
res(end+1,1) = isequal(I,I_);

% multiple rows and columns
c = [2 3 4; 5 6 7];
d = [1 0 1; 0 0 1];
intMat = intervalMatrix(c,d);
I = interval(intMat);
I_ = interval(c-d,c+d);
res(end+1,1) = isequal(I,I_);

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
