function P = minkDiff(P1,S,varargin)
% minkDiff - compute the Minkowski difference of two polytopes:
%         P1 - S = P <-> P + S \subseteq P1
%
% Syntax:  
%    P = minkDiff(P1,S)
%    P = minkDiff(P1,S,type)
%
% Inputs:
%    P1 - mptPolytope object
%    S - mptPolytope object, contSet object, or numerical vector
%    type - type of computation ('exact' or 'inner')
%
% Outputs:
%    P - mptPolytope object after Minkowski difference
%
% Example: 
%    P1 = mptPolytope([1 0 -1 0 1;0 1 0 -1 1]',[4;4;4;4;4]);
%    P2 = mptPolytope([-1 0 1;-2 3 0]');
%
%    P = minkDiff(P1,P2);
%
%    figure; hold on;
%    plot(P1);
%    plot(P2,[1,2],'r');
%    plot(P,[1,2],'g');
%
% References:
%    [1] M. Althoff, "On Computing the Minkowski Difference of Zonotopes"
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/minkDiff

% Author:       Niklas Kochdumper, Mark Wetzlinger
% Written:      04-February-2021
% Last update:  09-November-2022 (MW, rename 'minkDiff')
% Last revision:23-February-2023 (MW, update input argument handling,
%                                     different method for P-Z)

%------------- BEGIN CODE --------------

% default value
type = setDefaultValues({'exact'},varargin);
% check input arguments
inputArgsCheck({{P1,'att','mptPolytope'} ...
    {S,'att',{'contSet','numeric'}}...
    {type,'str',{'exact','inner'}}});

% different algorithms for different set representations
if isnumeric(S)
    P = P1 + (-S);

% exact computation 
elseif strcmp(type,'exact')
    
    if isa(S,'zonotope') || isa(S,'interval')

        % reduce offset of each halfspace of the polytope to obtain
        % inner-approximation of the Minkowski difference
        A = P1.P.A;
        b = P1.P.b;
        
        for i = 1:size(A,1)
            l = supportFunc_(S,A(i,:)','upper');
            b(i) = b(i) - l;
        end
        
        % instantiate polytope
        P = mptPolytope(A,b);
        
    elseif isa(S,'conZonotope') || isa(S,'mptPolytope') || ...
           isa(S,'zonoBundle')

        % compute Minkowski difference according to Lemma 1 in [1]
        V = vertices(S); P = P1 + (-V(:,1));

        for i = 2:size(V,2)
            P = and_(P,P1 + (-V(:,i)),'exact'); 
        end
        
    else
        throw(CORAerror('CORA:noExactAlg',P1,S));
    end

% inner-approximative computation 
elseif strcmp(type,'inner')
    % note: exact for all convex set representations, for non-convex an
    % inner-approximation?

    otherOptions = {};
    if isa(S,'polyZonotope') || isa(S,'conPolyZono')
        otherOptions = {'interval',8,1e-3};
    end
    
    % scale each halfspace of the polytope to obtain
    % inner-approximation of the Minkowski difference
    A = P1.P.A;
    b = P1.P.b;
    
    for i = 1:size(A,1)
        l = supportFunc_(S,A(i,:)','upper',otherOptions{:});
        b(i) = b(i) - l;
    end
    
    % instantiate polytope
    P = mptPolytope(A,b);          
    
end

%------------- END OF CODE --------------