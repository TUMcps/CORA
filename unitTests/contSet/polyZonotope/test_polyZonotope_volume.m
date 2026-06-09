function res = test_polyZonotope_volume
% test_polyZonotope_volume - unit test function for the volume of a
%   polynomial zonotope
%
% Syntax:
%    res = test_polyZonotope_volume
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
% See also: none

% Authors:       Niklas Kochdumper
% Written:       07-October-2025
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

    tol = 0.1;

    % number of factors = number of dimensions
    pZ = polyZonotope([0;0],[2 1 2;0 2 2],[],[1 0 2;0 1 1]);
    vol = volume(pZ);

    pgon = polygon(pZ,12);
    vol_ = volume(pgon);

    assert(vol_ > vol);
    assert(abs(vol - vol_) < tol);

    % test completed
    res = true;

% ------------------------------ END OF CODE ------------------------------
