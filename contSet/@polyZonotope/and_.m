function res = and_(pZ,S,varargin)
% and_ - overloads '&' operator, computes the intersection of a polynomial
%    zonotope and another set or numerical vector
%
% Syntax:
%    res = pZ & S
%    res = and_(pZ,S)
%
% Inputs:
%    pZ - polyZonotope object
%    S - contSet object
%
% Outputs:
%    res - intersection
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/and

% Authors:       Mark Wetzlinger
% Written:       28-September-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% call function with lower precedence
if isa(S,'contSet') && S.precedence < pZ.precedence
    res = and_(S,pZ,varargin{:});
    return
end

% is overridden in subclass if implemented; throw error
throw(CORAerror("CORA:noops",pZ,S));

% ------------------------------ END OF CODE ------------------------------
