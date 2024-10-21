function res = isIntersecting_(E,S,type,tol,varargin)
% isIntersecting_ - determines if an ellipsoid intersects a set
%
% Syntax:
%    res = isIntersecting_(E,S,type,tol)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object or double matrix (ncols = number of points)
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    E1 = ellipsoid([5 7;7 13],[1;2]);
%    E2 = ellipsoid(0.3*eye(2),[1;0]);
%    E3 = ellipsoid(0.3*eye(2),[2;0]);
%    Z = zonotope([3;0],0.5*eye(2));
%
%    isIntersecting(E1,E2)
%    isIntersecting(E1,E3)
%    isIntersecting(E1,Z,'approx')
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E2,[1,2],'g');
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(E3,[1,2],'r');
%
%    figure; hold on
%    plot(E1,[1,2],'b');
%    plot(Z,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, conZonotope/isIntersecting_

% Authors:       Victor Gassmann, Niklas Kochdumper
% Written:       13-March-2019 
% Last update:   21-November-2019 (NK, extended to other sets)
%                10-March-2021 (refactored, simplified)
%                19-May-2022 (VG, removed try-catch)
%                04-July-2022 (VG, adapted to class array definitions of some functions)
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[E,S] = reorderNumeric(E,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < E.precedence
    res = isIntersecting_(S,E,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(E,S,type,tol);
    return
end

% sets must not be empty
if representsa_(E,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% ellipsoid is just a point, check containment
if representsa_(E,'point',tol)
    res = contains_(S,E.q,type,tol);
    return;
end

% general method: use shortest distance (must be 0 for intersection)
try
    % we cannot check whether "distance(E,S)" is supported for some class S
    dist = distance(E,S);
    res = dist < E.TOL | withinTol(dist,E.TOL);
catch ME
    % distance not implemented for S
    if strcmp(type,'exact')
        throw(CORAerror('CORA:noExactAlg',E,S));
    end
    
    % use fallback: it is not certain that S implements quadMap
    if ismethod(S,'quadMap') && ismethod(E,'interval')
        res = priv_isIntersectingMixed(E,S);
    else
        throw(CORAerror('CORA:noops',E,S));
    end
end

% ------------------------------ END OF CODE ------------------------------
