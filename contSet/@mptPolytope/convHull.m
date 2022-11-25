function res = convHull(P1,varargin)
% convHull - computes the convex hull of two polytopes
%
% Syntax:  
%    res = convHull(P1,P2)
%
% Inputs:
%    P1 - first mptPolytope object
%    P2 - second mptPolytope object
%
% Outputs:
%    res - mptPolytope enclosing the convex hull of P1 and P2
%
% Example: 
%    poly1 = mptPolytope([-1 -1; 1 0;-1 0; 0 1; 0 -1],[2;3;2;3;2]);
%    poly2 = mptPolytope([-1 -1; -1 1; 1 1;0 -1],[0;2;2;0]) + [4;3];
%
%    poly = convHull(poly1,poly2);
%
%    figure
%    hold on
%    plot(poly1,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(poly2,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(poly,[1,2],'g');
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
       res = P1; return; 
    else
       P2 = varargin{1}; 
    end

    % find a polytope object
    if ~isa(P1,'mptPolytope')
        temp = P1;
        P1 = P2;
        P2 = temp;
    end

    % different cases depending on the class of the summand
    if isnumeric(P2)
        
        V2 = P2;

    elseif isa(P2,'mptPolytope') || isa(P2,'interval') || ...
           isa(P2,'conZonotope') || isa(P2,'zonoBundle') || ...
           isa(P2,'zonotope')

        % compute vertices
        V2 = vertices(P2);

    elseif isa(summand,'polyZonotope')

        res = P2 + P1;
        return;

    else
        % throw error for given arguments
        error(noops(P1,P2));
    end
    
    % comptue vertices of first polytope
    V1 = vertices(P1);
    
    % compute convex hull
    V = [V1,V2];
    
    K = convhulln(V');
    indices = unique(K);
    
    res = mptPolytope(V(:,indices)');

end

%------------- END OF CODE --------------