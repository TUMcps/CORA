function res = ne(pHA1,pHA2,varargin)
% ne - overloads '~=' operator to check if two parallel hybrid automata are
%    not equal
%
% Syntax:
%    res = pHA1 ~= pHA2
%    res = ne(pHA1,pHA2)
%    res = ne(pHA1,pHA2,tol)
%
% Inputs:
%    pHA1 - parallelHybridAutomaton object
%    pHA2 - parallelHybridAutomaton object
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

res = ~isequal(pHA1,pHA2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
