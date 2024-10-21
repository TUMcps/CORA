function I = interval(pgon)
% interval - compute an interval enclosure of the polygon
%
% Syntax:
%    I = interval(pgon)
%
% Inputs:
%    pgon - polygon
%
% Outputs:
%    I - interval
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

V = vertices_(pgon);
I = interval(min(V, [], 2), max(V, [], 2));

end

% ------------------------------ END OF CODE ------------------------------
