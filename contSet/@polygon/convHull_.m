function S_out = convHull_(pgon, S, varargin)
% convHull_ - compute convex hull of a polygon
%
% Syntax:
%    pgon = convHull_(pgon, S, varargin)
%
% Inputs:
%    pgon - polygon
%    S - contSet or numeric
%
% Outputs:
%    S_out - contSet    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/convHull

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, contSet integration)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

if nargin == 1
    S_out = polygon(convhull(pgon.set));
    return
end

% ensure that numeric is second input argument
[pgon,S] = reorderNumeric(pgon,S);

% check dimensions
equalDimCheck(pgon,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pgon.precedence
    S_out = convHull(S,pgon,varargin{:});
    return
end

% convert S into polygon
S_pgon = polygon(S);

% compute union
pgon_union = pgon | S_pgon;

% compute convex hull
S_out = polygon(convhull(pgon_union.set));

end

% ------------------------------ END OF CODE ------------------------------
