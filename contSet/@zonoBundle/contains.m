function res = contains(zB,S,varargin)
% contains - determines if a zonotope bundle contains a set or a point
%
% Syntax:  
%    res = contains(zB,S)
%    res = contains(zB,S,type)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object or single point
%    type - type of containment check ('exact' or 'approx')
%
% Outputs:
%    res - true/false
%
% Example: 
%    I1 = interval([0;-1],[2;2]);
%    I2 = I1 + [2;0];
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    contains(zB,I1)
%    contains(zB,I2)
%
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(I1,[1,2],'g');
%    
%    figure; hold on;
%    plot(zB,[1,2],'b');
%    plot(I2,[1,2],'r');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval/contains, conZonotope/contains

% Author:       Niklas Kochdumper
% Written:      19-November-2019
% Last update:  15-November-2022 (MW, return logical array for points)
%               25-November-2022 (MW, rename 'contains')
% Last revision:---

%------------- BEGIN CODE --------------

% pre-processing
[res,vars] = pre_contains('zonoBundle',zB,S,varargin{:});

% check premature exit
if res
    % if result has been found, it is stored in the first entry of var
    res = vars{1}; return
else
    zB = vars{1}; S = vars{2}; type = vars{3};
end


% point or point cloud in zonotope bundle containment
if isnumeric(S)
    
    res = false(1,size(S,2));
    for i = 1:size(S,2)
        res(i) = contains(conZonotope(zB),S(:,i));
    end
    
% capsule/ellipsoid in zonotope bundle containment
elseif isa(S,'capsule') || isa(S,'ellipsoid')
    
    P = mptPolytope(zB);
    res = contains(P,S); 

else
    
    % use the fast but over-approximative or the exact but possibly
    % slow containment check
    if strcmp(type,'exact')

        if isa(S,'taylm') || isa(S,'polyZonotope')
            throw(CORAerror('CORA:noExactAlg',S,"'taylm' or 'polyZonotope'"));
        elseif isa(S,'interval')
            res = contains(zB,vertices(S));
        else
            P = mptPolytope(zB);
            res = contains(P,S); 
        end
        
    else
        
        if isa(S,'taylm') || isa(S,'polyZonotope')
            P = mptPolytope(zB);
            res = contains(P,S); 
        else
            cZ = conZonotope(zB);
            res = contains(cZ,S,type);
        end
    end
end

%------------- END OF CODE --------------