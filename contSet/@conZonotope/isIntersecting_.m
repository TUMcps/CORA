function res = isIntersecting_(cZ,S,type,tol,varargin)
% isIntersecting_ - determines if a constrained zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(cZ,S)
%    res = isIntersecting_(cZ,S,type,tol)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - (optional) tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    % generate constrained zonotopes
%    Z = [0 2 -2 1;0 1.5 1 -1.5];
%    A = [1 1 1]; b = 1;
%    cZ1 = conZonotope(Z,A,b);
% 
%    Z = [1 2 0 0;1 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ2 = conZonotope(Z,A,b);
%
%    Z = [3 2 0 0;4 1 1 0];
%    A = [1 1 -1]; b = 0;
%    cZ3 = conZonotope(Z,A,b);
%
%    % check for intersection
%    isIntersecting(cZ1,cZ2)
%    isIntersecting(cZ1,cZ3)
%
%    % visualization
%    figure; hold on;
%    plot(cZ1,[1,2],'b');
%    plot(cZ2,[1,2],'g');
%
%    figure; hold on;
%    plot(cZ1,[1,2],'b');
%    plot(cZ3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, zonotope/isIntersecting_

% Authors:       Niklas Kochdumper
% Written:       21-November-2019
% Last update:   ---
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[cZ,S] = reorderNumeric(cZ,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < cZ.precedence
    res = isIntersecting_(S,cZ,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(cZ,S,type,tol);
    return
end

% % sets must not be empty (LPs too costly...)
% if representsa_(cZ,'emptySet',0) || representsa_(S,'emptySet',0)
%     res = false;
%     return
% end

% exact algorithm: check for non-empty intersection
if isa(S,'contSet') && strcmp(type,'exact')

    if isa(S,'conZonotope')
        res = ~representsa_(and_(cZ,S,'exact'),'emptySet',eps);
        return
    end
    if isa(S,'zonotope') || isa(S,'interval') || isa(S,'zonoBundle')
        res = ~representsa_(and_(cZ,conZonotope(S),'exact'),'emptySet',eps);
        return
    end
    
    throw(CORAerror('CORA:noExactAlg',cZ,S));
end
    
if isa(S,'contSet') && strcmp(type,'approx')
    if isa(S,'interval')
        res = isIntersecting_(polytope(S),cZ,type,tol);
    else
        res = isIntersecting_(polytope(interval(cZ)),S,type,tol);
    end
    return
end

throw(CORAerror('CORA:noops',cZ,S));

% ------------------------------ END OF CODE ------------------------------
