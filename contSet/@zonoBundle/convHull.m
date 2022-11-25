function res = convHull(zB1,varargin)
% convHull - computes the convex hull of two zonotope bundles
%
% Syntax:  
%    res = convHull(zB1,zB2)
%
% Inputs:
%    zB1 - first zonoBundle object
%    zB2 - second zonoBundle object
%
% Outputs:
%    res - zonoBundle enclosing the convex hull of zB1 and zB2
%
% Example: 
%    int = interval([4;3],[6;6]);
%    zono1 = zonotope([0 1 2 0;0 1 0 2]);
%    zono2 = zonotope([3 -0.5 3 0;-1 0.5 0 3]);
%    zB = zonoBundle({zono1,zono2});
%
%    res = convHull(int,zB);
%
%    figure
%    hold on
%    plot(int,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(zB,[1,2],'r','Filled',true,'EdgeColor','none');
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
        res = zB1; return;
    else
        zB2 = varargin{1}; 
    end

    % find a zonoBundle object
    if ~isa(zB1,'zonoBundle')
        temp = zB1;
        zB1 = zB2;
        zB2 = temp;
    end

    % different cases depending on the class of the second summand
    if isa(zB2,'zonoBundle') || isa(zB2,'interval') || ...
       isa(zB2,'zonotope') || isa(zB2,'conZonotope')

        poly2 = mptPolytope(zB2);

    elseif isa(zB2,'mptPolytope')

        poly2 = zB2;
        
    elseif isnumeric(zB2)
        
        poly2 = mptPolytope(zB2');

    elseif isa(summand,'polyZonotope') || isa(summand,'conPolyZono')

        res = zB2 + zB1;
        return

    else
        % throw error for given arguments
        error(noops(zB1,zB2));
    end
    
    % convert first zonoBundle to mptPolytope
    poly1 = mptPolytope(zB1);
    
    % compute convex hull
    poly = convHull(poly1,poly2);
    
    % convert to a zonotope bundle
    res = zonoBundle(poly);

end

%------------- END OF CODE --------------