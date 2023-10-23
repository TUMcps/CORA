function res = test_fullspace_lift
% test_fullspace_lift - unit test function of lift
%
% Syntax:
%    res = test_fullspace_lift
%
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
% Written:       06-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init fullspace
n = 4;
fs = fullspace(n);

% higher-dimensional space
N = 6;
dims = [1,2,5,6];
fs_ = lift(fs,N,dims);
fs_true = fullspace(N);
res = isequal(fs_true,fs_);

% not enough specified dimensions: fix error message before this...
% projDims = [1,2,5];
% try
%     fs_ = lift(fs,N,projDims);
%     res(end+1,1) = false;
% catch
%     res(end+1,1) = true;
% end

if CHECKS_ENABLED

% dimensions out of range
projDims = [-1,2,3,5];
try
    fs_ = lift(fs,N,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

% higher-dimensional space smaller than original space
N = 3;
dims = [1,2,5,6];
try
    fs_ = lift(fs,N,projDims);
    res(end+1,1) = false;
catch
    res(end+1,1) = true;
end

end

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
