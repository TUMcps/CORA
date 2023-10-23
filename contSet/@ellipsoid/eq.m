function res = eq(E1,E2,varargin)
% eq - Overloaded '==' operator for the comparison of ellipsoids
%
% Syntax:
%    res = E1 == E2
%    res = eq(E1,E2)
%    res = eq(E1,E2,tol)
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
%    res = E1==E2;
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ellipsoid/isequal

% Authors:       Victor Gassmann
% Written:       14-October-2019 
% Last update:   16-March-2021 (relative TOL)
%                04-July-2022 (VG, input checks)
%                23-December-2022 (MW, move code to isequal)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = isequal(E1,E2,varargin{:});

% ------------------------------ END OF CODE ------------------------------
