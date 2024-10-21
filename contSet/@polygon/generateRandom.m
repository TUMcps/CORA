function pgon = generateRandom(varargin)
% generateRandom - generate random polygon
%
% Syntax:
%    pgon = polygon.generateRandom()
%
% Inputs:
%    - 
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% generate random points
points = 10 * rand(1) * (-1 + 2 * rand(2, 100)) + 10 * (-1 + 2 * rand(2, 1));

% init polygon
pgon = polygon.enclosePoints(points);

end

% ------------------------------ END OF CODE ------------------------------
