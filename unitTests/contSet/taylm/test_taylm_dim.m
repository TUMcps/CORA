function res = test_taylm_dim
% test_taylm_dim - unit test function of dim
%
% Syntax:
%    res = test_taylm_dim
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

% Authors:       Tobias Ladner
% Written:       02-August-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

resvec = [];

% scalar case
tay = taylm(interval(1));
resvec(end+1) = dim(tay) == 1;
    
% higher-dimensional case
lb = [-3; -2; -5];
ub = [4; 2; 1];
I = interval(lb,ub);
tay = taylm(I);

% compute dimension
resvec(end+1) = dim(tay) == 3;

% combine results
res = all(resvec);

% ------------------------------ END OF CODE ------------------------------
