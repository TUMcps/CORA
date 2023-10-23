function res = eq(loc1,loc2,varargin)
% eq - overloads '==' operator to check if two locations objects are equal
%
% Syntax:
%    res = eq(loc1,loc2)
%    res = eq(loc1,loc2,tol)
%
% Inputs:
%    loc1 - location object
%    loc2 - location object
%    tol - tolerance (optional)
%
% Outputs:
%    res - true/false
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       10-January-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(loc1,loc2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
