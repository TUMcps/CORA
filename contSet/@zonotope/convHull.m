function res = convHull(Z1,varargin)
% convHull - computes an enclosure for the convex hull of two zonotopes
%
% Syntax:  
%    res = convHull(Z1,Z2)
%
% Inputs:
%    Z1 - first zonotope object
%    Z2 - second zonotope object
%
% Outputs:
%    res - zonotope enclosing the convex hull of Z1 and Z2
%
% Example: 
%    zono1 = zonotope([2 1 0; 2 0 1]);
%    zono2 = zonotope([-2 1 0; -2 0 1]);
%
%    zono = convHull(zono1,zono2);
%
%    figure; hold on;
%    plot(zono1,[1,2],'r','Filled',true,'EdgeColor','none');
%    plot(zono2,[1,2],'b','Filled',true,'EdgeColor','none');
%    plot(zono,[1,2],'g','LineWidth',3);
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
        res = Z1; return;
    else
        Z2 = varargin{1};
    end 

    % determine zonotope object
    if ~isa(Z1,'zonotope')
        temp = Z1;
        Z1 = Z2;
        Z2 = temp;
    end

    % different cases depending on the class of the second summand
    if isa(Z2,'zonotope') || isa(Z2,'interval') || isnumeric(Z2)

        % compute convex hull for constrained zonotopes
        cZ = convHull(conZonotope(Z1),conZonotope(Z2));
        
        % enclose result with a zonotope
        res = zonotope(cZ);

    elseif isa(summand,'mptPolytope') || isa(summand,'conZonotope') || ...
           isa(summand,'zonoBundle') || isa(summand,'polyZonotope') || ...
           isa(summand,'conPolyZono')

        res = Z2 + Z1;        

    else
        % throw error for given arguments
        error(noops(Z1,Z2));
    end
end

%------------- END OF CODE --------------