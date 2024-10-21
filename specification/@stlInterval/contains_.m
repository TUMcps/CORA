function res = contains_(I,S,type,tol,varargin)
% contains_ - determines if an STL interval contains a set or a point
%
% Syntax:
%    res = contains_(I,S)
%    res = contains_(I,S,type)
%    res = contains_(I,S,type,tol)
%
% Inputs:
%    I - interval object
%    S - contSet object or single point
%    type - 'exact' or 'approx'
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = stlInterval(1,5,false,false);
%    I2 = stlInterval(2,4);
%
%    contains(I1,I2)
%    contains(I1,5)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/contains

% Authors:       Florian Lercher
% Written:       06-February-2024
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

% set in empty set
if representsa_(I,'emptySet',0)
    res = representsa_(S,'emptySet',0);
    return
end

if isnumeric(S)
    % point in interval
    res = (I.lower < S && S < I.upper) || ...
        (I.leftClosed && withinTol(I.lower,S,tol)) || ...
        (I.rightClosed && withinTol(I.upper,S,tol));
elseif isa(S,'stlInterval')
    % interval in interval
    upperIncluded = I.upper > S.upper || ...
        (withinTol(I.upper,S.upper,tol) && (I.rightClosed || ~S.rightClosed));
    lowerIncluded = I.lower < S.lower || ...
        (withinTol(I.lower,S.lower,tol) && (I.leftClosed || ~S.leftClosed));
    res = lowerIncluded && upperIncluded;
elseif I.leftClosed && I.rightClosed
    res = contains(interval(I),S,type,tol);
elseif strcmp(type,'approx')
    if ~I.leftClosed
        lower = I.lower + tol + eps;
    else
        lower = I.lower;
    end
    if ~I.rightClosed
        upper = I.upper - tol - eps;
    else
        upper = I.upper;
    end
    res = contains(interval(lower,upper),S,type,tol);
else
    throw(CORAerror('CORA:noExactAlg',I,S));
end
    

% ------------------------------ END OF CODE ------------------------------
