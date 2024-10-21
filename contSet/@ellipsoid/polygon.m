function pgon = polygon(E,varargin)
% polygon - convert a set to a (inner-approximative) polygon
%
% Syntax:
%    pgon = polygon(E)
%
% Inputs:
%    E - contSet object
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
% Written:       14-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------
 
% check dimension of set
if dim(E) ~= 2
    throw(CORAerror('CORA:wrongValue','first','Given set must be 2-dimensional.'));
end

% compute boundary points
if representsa_(E,'point',eps)  % only center remaining
    % read center point in desired dimensions
    V = E.q;

else
    % number of points on boundary
    N = 1000;

    % compute points
    V = priv_boundary(E,N);
end

% init polygon
pgon = polygon(V);

end

% ------------------------------ END OF CODE ------------------------------
 