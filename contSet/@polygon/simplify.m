function pgon = simplify(pgon,varargin)
% simplify - enclose the polygon by a simpler polygon with less vertices
%
% Syntax:
%    pgon = simplify(pgon)
%    pgon = simplify(pgon,tol)
%
% Inputs:
%    pgon - polygon
%    tol - numeric, tolerance
%
% Outputs:
%    pgon - polygon
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Niklas Kochdumper
% Written:       13-March-2020
% Last update:   11-October-2024 (TL, moved to polygon/compact)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

CORAwarning('CORA:deprecated','function','polygon/simplify','CORA v2025', ...
    'When updating the code, please replace every function call ''simplify(polygon,tol)'' with ''compact(polygon,''douglasPeucker'',tol)''.', ...
    'This change was made in an effort to unify the syntax across all set representations.')

% set default values
tol = setDefaultValues({0.01},varargin);

% call compact
pgon = compact_(pgon,'douglasPeucker',tol);

end

% ------------------------------ END OF CODE ------------------------------
