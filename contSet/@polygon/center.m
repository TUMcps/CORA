function c = center(pgon)
% center - get the center of the polygon
%
% Syntax:
%    c = center(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    c - numeric
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

[x, y] = centroid(pgon.set);
c = [x; y];

end

% ------------------------------ END OF CODE ------------------------------
