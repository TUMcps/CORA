function res = test_emptySet_project
% test_emptySet_project - unit test function of project
%
% Syntax:
%    res = test_emptySet_project
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
% Written:       05-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init empty set
n = 4;
O = emptySet(n);

% project to subspace
projDims = [1,3];
O_ = project(O,projDims);

% true solution
O_true = emptySet(length(projDims));

% compare solutions
res = O_ == O_true;

% subspace out of range
projDims = [-1,2];
try
    O_ = project(O,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% subspace out of range
projDims = [3,5];
try
    O_ = project(O,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
