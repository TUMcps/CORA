function cZ = cartesianProduct(cZ,Z)
% cartesianProduct - Returns the cartesian product of a constrained 
%                    zonotope and a zonotope object
%
% Syntax:  
%    cZ = cartesianProduct(cZ,Z)
%
% Inputs:
%    cZ - conZonotope object
%    Z - zonotope object
%
% Outputs:
%    cZ - resulting conZonotope object
%
% Example: 
%    Z = [0, 1 2];
%    A = [1 1];
%    b = 1;
%    cZ = conZonotope(Z,A,b);
%    zono = zonotope([0 1]);
%
%    cZcart = cartesianProduct(cZ,zono);
%
%    plotFilled(cZcart,[1,2],'r');
%    xlim([0.5 2.5]);
%    ylim([-1.5 1.5]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: zonotope/cartesianProduct

% Author:       Niklas Kochdumper
% Written:      10-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

    % new center vector
    c = [cZ.Z(:,1); Z.Z(:,1)];

    % determine sizes of generator matrices
    [r1, c1] = size(cZ.Z(:,2:end));
    [r2, c2] = size(Z.Z(:,2:end));

    % new generator matrix
    G = [cZ.Z(:,2:end), zeros(r1,c2); zeros(r2,c1), Z.Z(:,2:end)];

    % new constraint matrix
    A = [cZ.A, zeros(size(cZ.A,1),r2)];

    % generate resulting constrained zonotope
    cZ.Z = [c,G];
    cZ.A = A;

    cZ.ksi = [];
    cZ.R = [];

%------------- END OF CODE --------------