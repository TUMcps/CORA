function E = enclose(varargin)
% enclose - Generates an ellipsoid that encloses an ellipsoid and its linear 
%           transformation according to Prop. 2.7 in [1]
%
% Syntax:  
%    E = enclose(E1,E2)
%    E = enclose(E1,M,Eplus)
%
% Inputs:
%    E1 - first ellipsoid object
%    E2 - second ellipsoid object, satisfying E2 = (M * E1) + Eplus
%    M - matrix for the linear transformation
%    Eplus - ellipsoid object added to the linear transformation
%
% Outputs:
%    E - ellipsoid that encloses E1 and E2
%
% Example: 
%    E1=ellipsoid(eye(2));
%    E2=ellipsoid([1,0;0,3],[1;-1]);
%    E=enclose(E1,E2);
%    plot(E1);
%    hold on
%    plot(E2);
%    plot(E);
%
% References:
%    [1] E. Yildirim. "On the minimum volume covering ellipsoid of
%        ellipsoids", 2006
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/enclose

% Author:       Victor Gassmann
% Written:      13-March-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % parse input arguments
    if nargin == 2
        E1 = varargin{1};
        E2 = varargin{2};
    else
        E1 = varargin{1};
        M = varargin{2};
        Eplus = varargin{3};

        E2 = (M*E1) + Eplus;
    end

    % compute enclosure using convex hull
    E = convHull(E1,E2);

%------------- END OF CODE --------------