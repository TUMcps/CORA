function res = isIntersecting_(I,S,type,tol,varargin)
% isIntersecting_ - determines if an interval intersects a set
%
% Syntax:
%    res = isIntersecting_(I,S,type,tol)
%
% Inputs:
%    I - interval object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
%    tol - tolerance
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([0;0],[2;2]);
%    I2 = interval([1;1],[3;3]);
%    I3 = interval([-3;-3],[-1;1]);
%
%    isIntersecting(I1,I2)
%    isIntersecting(I1,I3)
%
%    figure; hold on;
%    plot(I1,[1,2],'b');
%    plot(I2,[1,2],'g');
%
%    figure; hold on;
%    plot(I1,[1,2],'b');
%    plot(I3,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: contSet/isIntersecting, zonotope/isIntersecting_

% Authors:       Matthias Althoff, Niklas Kochdumper
% Written:       22-July-2016
% Last update:   14-September-2019
%                21-November-2019 (NK, added intersection with other sets)
%                12-March-2021 (MW, add empty case)
%                04-December-2023 (MW, fix empty case)
%                24-May-2024 (TL, added interval-numeric case)
% Last revision: 27-March-2023 (MW, rename isIntersecting_)

% ------------------------------ BEGIN CODE -------------------------------

% ensure that numeric is second input argument
[I,S] = reorderNumeric(I,S);

% call function with lower precedence
if isa(S,'contSet') && S.precedence < I.precedence
    res = isIntersecting_(S,I,type,tol);
    return
end

% numeric case: check containment
if isnumeric(S)
    res = contains_(I,S,type,tol);
    return
end

% sets must not be empty
if representsa_(I,'emptySet',0) || representsa_(S,'emptySet',0)
    res = false;
    return
end

% interval and interval intersection
if isa(S,'interval')
    
    % get object properties
    sup1 = I.sup; inf1 = I.inf;
    sup2 = S.sup; inf2 = S.inf;
    
    % loop over all dimensions separately
    for i = 1:length(I)
        if ~aux_isIntersecting1D(inf1(i),sup1(i),inf2(i),sup2(i),tol)
            res = false;
            return
        end
    end
    res = true;
    return
end

throw(CORAerror('CORA:noops',I,S));

% if strcmp(type,'approx')
%     res = isIntersecting_(polytope(I),S,type,tol);
% end
    
end


% Auxiliary functions -----------------------------------------------------

function res = aux_isIntersecting1D(inf1,sup1,inf2,sup2,tol)
% check if two one-dimensional intervals intersect
    res = false;

    if inf1 <= inf2 || withinTol(inf1,inf2,tol)
        if inf2 <= sup1 || withinTol(inf2,sup1,tol)
            res = true;
        end
        
    else % inf2 < inf1
        if inf1 <= sup2 || withinTol(inf1,sup2,tol)
            res = true;
        end
    end

end

% ------------------------------ END OF CODE ------------------------------
