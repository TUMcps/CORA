function pgon = subtract(pgon1, pgon2, varargin)
% subtract - subtracts the second polygon from the first
%
% Syntax:
%    pgon = subtract(pgon1, pgon2)
%
% Inputs:
%    pgon1 - polygon
%    pgon2 - polygon
%    varargin - additional parameters for polyshape/subtract
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: polyshape/subtract

% Authors:       Tobias Ladner
% Written:       13-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input
narginchk(2,inf);
inputArgsCheck({ ...
    {pgon1,'att','polygon'}; ...
    {pgon2,'att','polygon'}; ...
})

% call polyshape/subtract
pshape = subtract(pgon1.set,pgon2.set);

% convert back to polygon
pgon = polygon(pshape);

end

% ------------------------------ END OF CODE ------------------------------
