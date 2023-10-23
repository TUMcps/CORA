function zB = enclose(zB,varargin)
% enclose - encloses a zonotope bundle and its affine transformation (see
%    Proposition 5 in [1])
%
% Description:
%    Computes the set
%    { a x1 + (1 - a) * x2 | x1 \in zB, x2 \in S, a \in [0,1] }
%    where S = M*zB + Splus
%
% Syntax:
%    zB = enclose(zB,S)
%    zB = enclose(zB,M,Splus)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%    M - matrix for the linear transformation
%    Splus - zonoBundle object added to the linear transformation
%
% Outputs:
%    zB - zonoBundle object that encloses given zonotope bundle and set
%
% References:
%    [1] M. Althoff. "Zonotope bundles for the efficient computation of 
%        reachable sets", 2011
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Authors:       Matthias Althoff
% Written:       10-November-2010 
% Last update:   25-January-2016
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% parse input arguments
if nargin == 2
    S = varargin{1};
elseif nargin == 3
    % evaluate M*zB + Splus
    S = (varargin{1}*zB) + varargin{2};
else
    throw(CORAerror('CORA:tooManyInputArgs',3));
end


if isa(S,'zonoBundle')

    % compute enclosure for each zonotope pair
    for i=1:zB.parallelSets
        zB.Z{i} = enclose(zB.Z{i},S.Z{i});
    end

elseif isa(S,'zonotope')

    % compute enclosure for each zonotope
    for i=1:zB.parallelSets
        zB.Z{i} = enclose(zB.Z{i},S);
    end

else

    % not implemented for other contSet classes
    throw(CORAerror('CORA:noops',zB,S));

end

% ------------------------------ END OF CODE ------------------------------
