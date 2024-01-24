function res = isequal(Z1,Z2,varargin)
% isequal - checks if two zonotopes are equal (note: no deletion of aligned
%    generators since this is quite costly)
%
% Syntax:
%    res = isequal(Z1,Z2)
%    res = isequal(Z1,Z2,tol)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    Z1 = zonotope(zeros(2,1),[1 0 2; 2 0 1]);
%    Z2 = zonotope(zeros(2,1),[1 2 0; 2 1 0]);
%    isequal(Z1,Z2)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Authors:       Mark Wetzlinger
% Written:       16-September-2019
% Last update:   09-June-2020 (include re-ordering of generators)
%                13-November-2022 (MW, integrate modular check)
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% too many input arguments
if nargin > 3
    throw(CORAerror('CORA:tooManyInputArgs',3));
end

% parse input arguments
tol = setDefaultValues({eps},varargin);

% check input arguments
inputArgsCheck({{Z1,'att','zonotope'};
                {Z2,'att',{'zonotope','interval'}};
                {tol,'att','numeric',{'nonnan','scalar','nonnegative'}}});

% init result
res = false;

% convert interval to zonotope
if isa(Z2,'interval')
    Z2 = zonotope(Z2);
end

% compare dimensions (quick check)
if dim(Z1) ~= dim(Z2)
    return
end

% compare centers (quick check)
if ~all(withinTol(Z1.c,Z2.c,tol))
    return
end

% delete zeros from generator matrices
G1 = compact_(Z1,'zeros',tol).G;
G2 = compact_(Z2,'zeros',tol).G;

% compare number of generators
if size(G1,2) ~= size(G2,2)
    return
end

% compare generator matrices: must match with no remainder, order is
% irrelevant, sign may be inverted
res = compareMatrices(G1,G2,tol,'equal',false,false);

% ------------------------------ END OF CODE ------------------------------
