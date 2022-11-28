function zB = convHull(zB,varargin)
% convHull - computes the convex hull of two zonotope bundles
%
% Syntax:  
%    res = convHull(zB,S)
%
% Inputs:
%    zB - zonoBundle object
%    S - contSet object
%
% Outputs:
%    zB - zonoBundle enclosing the convex hull
%
% Example: 
%    I = interval([4;3],[6;6]);
%    Z1 = zonotope([0 1 2 0;0 1 0 2]);
%    Z2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({Z1,Z2});
%
%    res = convHull(zB,I);
%
%    figure; hold on;
%    plot(int,[1,2],'FaceColor','b');
%    plot(zB,[1,2],'FaceColor','r');
%    plot(res,[1,2],'g','LineWidth',3);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: conZonotope/convHull

% Author:        Niklas Kochdumper
% Written:       26-November-2019 
% Last update:   05-May-2020 (MW, standardized error message)
% Last revision: ---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin == 1
        return;
    else
        S = varargin{1}; 
    end

    % find a zonoBundle object
    if ~isa(zB,'zonoBundle')
        temp = zB;
        zB = S;
        S = temp;
    end

    % different cases depending on the class of the second summand
    if isa(S,'zonoBundle') || isa(S,'interval') || ...
        isa(S,'zonotope') || isa(S,'conZonotope')

        P2 = mptPolytope(S);

    elseif isa(S,'mptPolytope')

        P2 = S;
        
    elseif isnumeric(S)
        
        P2 = mptPolytope(S');

    elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')

        zB = S + zB; return;

    else
        % throw error for given arguments
        throw(CORAerror('CORA:noops',zB,S));
    end
    
    % convert first zonoBundle to mptPolytope
    P1 = mptPolytope(zB);
    
    % compute convex hull
    P = convHull(P1,P2);
    
    % convert to a zonotope bundle
    zB = zonoBundle(P);

end

%------------- END OF CODE --------------