function res = test_fullspace_project
% test_fullspace_project - unit test function of project
%
% Syntax:
%    res = test_fullspace_project
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

% init fullspace
n = 4;
fs = fullspace(n);

% project to subspace
projDims = [1,4];
fs_ = project(fs,projDims);
fs_true = fullspace(length(projDims));
res = isequal(fs_true,fs_);

% subspace out of range
projDims = [-1,2];
try
    fs_ = project(fs,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% subspace out of range
projDims = [3,5];
try
    fs_ = project(fs,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
