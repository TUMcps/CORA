function res = isequal(pgon1, pgon2, varargin)
% isequal - compares whether two polygons represent the same shape
% notes: - (almost) collinear points are removed by the constructor
%        - vertices have to be in order, but start of list can vary
%
% Syntax:
%    res - isequal(pgon1,pgon2,varargin)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%    tol - numeric, tolerance
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       27-May-2023
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set tolerance
tol = setDefaultValues({1e-8}, varargin);

% the constructor already removes collinear points, but within its
% own tolerance -> look for collinear points up to given tolerance
pgon1 = removeCollinearPoints(pgon1, tol);
pgon2 = removeCollinearPoints(pgon2, tol);

% read out vertices
V1 = vertices_(pgon1)';
V2 = vertices_(pgon2)';

% check empty case
if isempty(V1) && isempty(V2)
    res = true; return
end

% since the constructor removes collinear points, the number of
% vertices has to be equal
if size(V1, 1) ~= size(V2, 1)
    res = false; return
end

% to deal with the different start vertex, first find one vertex
% that is part of both lists
% -> check all occurrences of the first vertex of P1 in P2
idx = find(all(withinTol(V1(1, :), V2, tol), 2));

% any matching vertex found?
if isempty(idx)
    res = false; return
end

% loop over all potential start vertices
for i = 1:length(idx)
    % re-order vertices of P2
    V2_ = [V2(idx:end, :); V2(1:idx-1, :)];
    % compare in order
    if compareMatrices(V1, V2_, tol, 'equal', true)
        res = true; return
    end
end

% no re-ordering successful
res = false;

end

% ------------------------------ END OF CODE ------------------------------
