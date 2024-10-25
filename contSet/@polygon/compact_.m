function pgon = compact_(pgon,method,tol)
% compact_ - removes redundancies in the representation of a set
%
% Syntax:
%    pgon = compact_(pgon)
%    pgon = compact_(pgon, method, tol)
%
% Inputs:
%    pgon - polygon object
%    method - method
%    tol - tolerance
%
% Outputs:
%    pgon - polygon object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/compact

% Authors:       Tobias Ladner
% Written:       11-October-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% enlarge pgon by tolerance
pgon = expandBoundaries(pgon, tol);

% apply further simplifications based on chosen method
switch method
    case 'douglasPeucker'
        pgon = aux_douglasPeucker(pgon,tol);
    case 'simplify'
        pgon = aux_simplify(pgon,tol);
    case 'all'
        pgon = aux_douglasPeucker(pgon,tol);
        pgon = aux_simplify(pgon,tol);
end

end


% Auxiliary functions -----------------------------------------------------

function pgon = aux_douglasPeucker(pgon,tol)
    % simplify polygon boundary using the Douglas-Peucker algorithm
    V = vertices_(Vertices);
    V_ = douglasPeucker(V', tol);
    pgon = polygon(V_(1, :), V_(2, :));
end

function pgon = aux_simplify(pgon,tol)
    pgon = polygon(simplify(pgon.set));
end

% ------------------------------ END OF CODE ------------------------------
