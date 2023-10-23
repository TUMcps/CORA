function res = eq(HA1,HA2,varargin)
% eq - overloads '==' operator to check if two hybrid automata are equal
%
% Syntax:
%    res = eq(HA1,HA2)
%    res = eq(HA1,HA2,tol)
%
% Inputs:
%    HA1 - hybridAutomaton object
%    HA2 - hybridAutomaton object
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

res = isequal(HA1,HA2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
