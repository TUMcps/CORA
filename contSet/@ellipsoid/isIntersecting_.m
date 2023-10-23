function res = isIntersecting_(E,S,type,varargin)
% isIntersecting_ - determines if an ellipsoid intersects a set
%
% Syntax:
%    res = isIntersecting_(E,S)
%    res = isIntersecting_(E,S,type)
%
% Inputs:
%    E - ellipsoid object
%    S - contSet object or double matrix (ncols = number of points)
%    type - type of check ('exact' or 'approx')
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

% all-zero shape matrix -> ellipsoid is just a point, check containment
if rank(E)==0
    res = contains_(S,E.q,type,100*eps);
    return;
end

if isa(S,'double')
    res = contains_(E,S,type,100*eps);
    return;
else
    % make sure that contSet object is scalar
    inputArgsCheck({{S,'att',{'contSet'},{'scalar'}}});
end

try
    % necessary since while we can check if "distance" exists for E, we  
    % cannot check whether "distance(E,S)" for S of some class is supported
    % use distance
    tmp = distance(E,S);
    res = tmp < E.TOL | withinTol(tmp,E.TOL);
catch 
    % distance not implemented for S
    % use fallback
    % exact or over-approximative algorithm
    if strcmp(type,'exact')
        throw(CORAerror('CORA:noops',E,S));
    end
    
    % it is not certain that S implements quadMap
    if ismethod(S,'quadMap') && ismethod(E,'interval')
        res = isIntersectingMixed(E,S);
    else
        throw(CORAerror('CORA:noops',E,S));
    end
end

% ------------------------------ END OF CODE ------------------------------
