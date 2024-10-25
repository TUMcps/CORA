function pgon = polygon(S,varargin)
% polygon - convert a set to a (outer-approximative) polygon
%
% Syntax:
%    pgon = polygon(S)
%
% Inputs:
%    S - contSet object
%    varargin - additonal parameters for vertices computation
%
% Outputs:
%    pgon - polygon object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polygon, contSet/vertices

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% check dimension of set
if dim(S) ~= 2
    throw(CORAerror('CORA:wrongValue','first','Given set must be 2-dimensional.'));
end

% compute vertices
V = vertices(S,varargin{:});

% init polygon
pgon = polygon(V);

try
    if representsa(S,'convexSet')
        pgon = convHull_(pgon);
    end
catch ME
    % not always implemented / hard to determine for some sets
    % keep pgon as is
end

end

% ------------------------------ END OF CODE ------------------------------
 