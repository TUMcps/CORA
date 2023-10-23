function res = isIntersecting_(cZ,S,type,varargin)
% isIntersecting_ - determines if a constrained zonotope intersects a set
%
% Syntax:
%    res = isIntersecting_(cZ,S)
%    res = isIntersecting_(cZ,S,type)
%
% Inputs:
%    cZ - conZonotope object
%    S - contSet object
%    type - type of check ('exact' or 'approx')
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

% intersection with halfspace, hyperplane or polytope
if isa(S,'halfspace') || isa(S,'conHyperplane') || ...
   isa(S,'polytope') || isa(S,'ellipsoid')

    res = isIntersecting_(S,cZ,type);

else
    
    % exact or over-approximative algorithm
    if strcmp(type,'exact')
        
        % convert objects to constrained zonotopes
        if isa(S,'zonotope') || isa(S,'interval')
           S = conZonotope(S); 
        end
        
        % conZonotope and conZonotope intersection
        if isa(S,'conZonotope')
           res = ~representsa_(and_(cZ,S,'exact'),'emptySet',eps);
        elseif isa(S,'zonoBundle')
           res = isIntersecting_(S,cZ,type); 
        else
            throw(CORAerror('CORA:noops',cZ,S));
        end
        
    else       
        if isa(S,'inverval')
           res = isIntersecting_(S,cZ,type); 
        else
           res = isIntersecting_(interval(cZ),S,type); 
        end
    end
end

% ------------------------------ END OF CODE ------------------------------
