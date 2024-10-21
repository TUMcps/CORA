function han = plot3D(C,NVpairsPlot,varargin)
% plot3D - plots a 3D projection of a capsule
%
% Syntax:
%    han = plot3D(C)
%    han = plot3D(C,dims)
%    han = plot3D(C,dims,NVpairs)
%
% Inputs:
%    C - capsule object
%    dims - (optional) dimensions for projection
%    NVpairsPlot - (optional) plot settings (LineSpec and Name-Value pairs)
%
% Outputs:
%    han - handle to the graphics object
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot

% Authors:       Matthias Althoff
% Written:       04-March-2019
% Last update:   ---
% Last revision: 14-October-2024 (TL, split plot1D/plot2D/plot3D)

% ------------------------------ BEGIN CODE -------------------------------

if isempty(C.g)
    % generator is empty

    E = ellipsoid(eye(3)*C.r^2, C.c);
    Z = zonotope(E, 100, 'outer:norm');
    han = plot(Z, 1:3, NVpairsPlot{:});

else
    % generator is not empty

    % generate ellipsoids for spheres
    E1 = ellipsoid(eye(3)*C.r^2, C.c+C.g);
    E2 = ellipsoid(eye(3)*C.r^2, C.c-C.g);

    % enclose ellipsoid with a zonotope
    Z1 = zonotope(E1,'outer:norm',100);
    Z2 = zonotope(E2,'outer:norm',100);

    % compute vertices of the zonotopes
    V1 = vertices(Z1);
    V2 = vertices(Z2);

    % plot the convex hull
    han = plotPolytope3D([V1, V2], NVpairsPlot{:});
end

end

% ------------------------------ END OF CODE ------------------------------
