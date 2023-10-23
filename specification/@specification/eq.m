function res = eq(spec1,spec2,varargin)
% eq - overloads '==' operator for comparison of specification objects
%
% Syntax:
%    res = spec1 == spec2
%    res = eq(spec1,spec2)
%    res = eq(spec1,spec2,tol)
%
% Inputs:
%    spec1 - specification object
%    spec2 - specification object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       29-April-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(spec1,spec2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
