function res = ne(C1,C2,varargin)
% ne - overloaded '~=' operator for exact comparison of two capsules
%
% Syntax:
%    res = C1 ~= C2
%    res = ne(C1,C2)
%    res = ne(C1,C2,tol)
%
% Inputs:
%    C1 - capsule object
%    C2 - capsule object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    C1 = capsule([1; 1; 0], [0.5; -1; 1], 0.5);
%    C2 = capsule([1; 0; 0], [0.5; -1; 1], 0.5);
%    C1 ~= C2
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: capsule/isequal

% Authors:       Mark Wetzlinger
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(C1,C2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
