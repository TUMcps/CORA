function res = contains_(pgon1, pgon2, method, tol, varargin)
% contains_ - check if polygon object pgon1 contains polygon pgon2 or a vector 
%    of points (one logical value for each point)
%
% Syntax:
%    res = contains_(pgon1, pgon2, varargin)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%    method - 'exact'
%    tol - numeric, tolerance
%
% Outputs:
%    res - logical
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, integration in contSet)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enlarge pgon1 by tolerance
pgon1 = expandBoundaries(pgon1, tol);

if isnumeric(pgon2)

    res = isinterior(pgon1.set, pgon2(1, :), pgon2(2, :))';

elseif isa(pgon2, 'contSet')
    % convert to polygon
    pgon2 = polygon(pgon2);

    % compute union
    u = union(pgon1.set, pgon2.set);

    % check if area of pgon1 is identical to area of union
    A1 = area(pgon1.set);
    A2 = area(u);

    res = withinTol(A1, A2, tol);
else
    throw(CORAerror('CORA:notSupported', ...
        'This set representation is not supported!'));
end

end

% ------------------------------ END OF CODE ------------------------------
