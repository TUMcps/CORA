function [res, S] = representsa_(pgon, type, tol, varargin)
% representsa_ - check if the given polygon represents another set
%
% Syntax:
%    res = representsa_(pgon,type)
%    [res,S] = representsa_(pgon,type)
%
% Inputs:
%    pgon - polygon
%    type - char, type which should be checked
%    tol - tolerance
%
% Outputs:
%    res - logical
%    S - converted set
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/representsa

% Authors:       Mark Wetzlinger, Tobias Ladner
% Written:       13-March-2020
% Last update:   16-October-2024
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% init second output argument (covering all cases with res = false)
S = [];

% check type
switch type
    case 'emptySet'
        res = isempty(vertices_(pgon));
        S = emptySet(2);

    case 'fullspace'
        res = false;

    case 'polygon'
        res = true;
        S = pgon;

    case 'point'
        V = vertices_(pgon);
        res = size(V,2) == 1;
        if res
            S = V(:,1);
        end

    case 'origin'
        res = all(withinTol(vertices_(pgon),0,tol),"all");
        S = [0;0];

    otherwise
        throw(CORAerror('CORA:notSupported', sprintf('Type ''%s'' is not supported.',type)));

end

end

% ------------------------------ END OF CODE ------------------------------
