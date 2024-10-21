function pgon = or(pgon1, pgon2, varargin)
% or - computes the union of two polygons
%
% Syntax:
%    pgon = or(pgon1, pgon2, varargin)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%    varargin - additional parameters for union of polyshape
%
% Outputs:
%    han - plot handle
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


if representsa_(pgon1, 'emptySet', eps)
    pgon = pgon2;
elseif representsa_(pgon2, 'emptySet', eps)
    pgon = pgon1;
else
    pgon = polygon(union(pgon1.set, pgon2.set, varargin{:}));
end

end

% ------------------------------ END OF CODE ------------------------------
