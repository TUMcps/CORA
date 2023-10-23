function res = eq(fs,S,varargin)
% eq - overloads '=='-operator for full-dimensional spaces
%    case R^0: see isequal
%
% Syntax:
%    res = fs == S
%    res = eq(fs,S)
%    res = eq(fs,S,tol)
%
% Inputs:
%    fs - fullspace object
%    S - contSet object, numerical vector
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    fs1 = fullspace(2);
%    fs2 = fullspace(3);
%    res1 = fs1 == fs1;
%    res2 = fs1 == fs2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       22-March-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% note: tolerance has no effect, only for overloading purposes
res = isequal(fs,S,varargin{:});

% ------------------------------ END OF CODE ------------------------------
