function res = isIntersecting(cZ,S,varargin)
% isIntersecting - determines if a constrained zonotope intersects a set
%
% Syntax:  
%    res = isIntersecting(cZ,S)
%    res = isIntersecting(cZ,S,type)
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
% See also: zonotope/isIntersecting

% Author:       Niklas Kochdumper
% Written:      21-Nov-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[resFound,vars] = pre_isIntersecting('conZonotope',cZ,S,varargin{:});

% check premature exit
if resFound
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    % assign values
    cZ = vars{1}; S = vars{2}; type = vars{3};
end


% intersection with halfspace, hyperplane or polytope
if isa(S,'halfspace') || isa(S,'conHyperplane') || ...
   isa(S,'mptPolytope') || isa(S,'ellipsoid')

    res = isIntersecting(S,cZ,type);

else
    
    % exact or over-approximative algorithm
    if strcmp(type,'exact')
        
        % convert objects to constrained zonotopes
        if isa(S,'zonotope') || isa(S,'interval')
           S = conZonotope(S); 
        end
        
        % conZonotope and conZonotope intersection
        if isa(S,'conZonotope')
           res = ~isempty(cZ & S);
        elseif isa(S,'zonoBundle')
           res = isIntersecting(S,cZ,type); 
        else
            throw(CORAerror('CORA:noops',cZ,S));
        end
        
    else       
        if isa(S,'inverval')
           res = isIntersecting(S,cZ,type); 
        else
           res = isIntersecting(interval(cZ),S,type); 
        end
    end
end

%------------- END OF CODE --------------