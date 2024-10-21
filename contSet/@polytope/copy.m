function P_out = copy(P)
% copy - copies the polytope object (used for dynamic dispatch)
%
% Syntax:
%    P_out = copy(P)
%
% Inputs:
%    P - polytope object
%
% Outputs:
%    P_out - copied polytope object
%
% Example: 
%    P = polytope([1 0; -1 1; -1 -1],[1;1;1]);
%    P_out = copy(P);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       30-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call copy constructor
P_out = polytope(P);

% ------------------------------ END OF CODE ------------------------------
