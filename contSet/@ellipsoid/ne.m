function res = ne(E1,E2,varargin)
% ne - Overloaded '~=' operator for the comparison of ellipsoids
%
% Syntax:
%    res = E1 ~= E2
%    res = ne(E1,E2)
%    res = ne(E1,E2,tol)
%
% Inputs:
%    E1 - ellipsoid object 
%    E2 - ellipsoid object 
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([1 0; 0 1]);
%    E2 = ellipsoid(E1.Q,E1.q);
%    res = E1 ~= E2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/isequal

% Authors:       Mark Wetzlinger
% Written:       23-December-2022
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = ~isequal(E1,E2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
