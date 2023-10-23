function res = test_fullspace_box
% test_fullspace_box - unit test function of box
%
% Syntax:
%    res = test_fullspace_box
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
n = 2;
fs = fullspace(n);

% compute box
fs_ = box(fs);

% compare results
res = isequal(fs,fs_);

% ------------------------------ END OF CODE ------------------------------
