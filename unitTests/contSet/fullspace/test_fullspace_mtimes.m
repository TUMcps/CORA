function res = test_fullspace_mtimes
% test_fullspace_mtimes - unit test function of mtimes
%
% Syntax:
%    res = test_fullspace_mtimes
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
n = 3;
fs = fullspace(n);

% multiplication with scalar
p = 2;
fs_ = fs * p;
res = isequal(fs,fs_);
fs_ = p * fs;
res(end+1,1) = isequal(fs,fs_);

% multiplication with full-rank square matrix (still fullspace)
M = [2 1 0; -1 -1 2; 1 2 1];
fs_ = M * fs;
res(end+1,1) = isequal(fs,fs_);

% multiplication as projection
M = [2 0 0; 0 1 0; 0 0 0];
fs_ = M * fs;
res(end+1,1) = fs_ == interval([-Inf;-Inf;0],[Inf;Inf;0]);

% multiplication as projection onto higher-dimensional space
M = [eye(n); [1 0 0]];
fs_ = M * fs;
res(end+1,1) = isequal(fs_,fullspace(n+1));

% combine results
res = all(res);

% ------------------------------ END OF CODE ------------------------------
