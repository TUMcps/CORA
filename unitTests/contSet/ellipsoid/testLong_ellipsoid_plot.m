function res = testLong_ellipsoid_plot
% testLong_ellipsoid_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check
%
% Syntax:
%    res = testLong_ellipsoid_plot
%
% Inputs:
%    -
%
% Outputs:
%    res - true/false
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Authors:       Mark Wetzlinger
% Written:       04-August-2020
% Last update:   ---
% Last revision: ---

% ------------------------------ BEGIN CODE -------------------------------

res = true;

% create a random psd Q matrix
dim_x = 4;
O = orth(randn(dim_x));
D = diag(abs(randn(dim_x,1)) + 0.3);
Q = O*D*O';

% generate ellipsoid
E = ellipsoid(Q);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(E);
    
    % two arguments: object, dimensions
    plot(E,[1,2]);
    plot(E,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(E,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(E,[1,2],'LineWidth',2);
    plot(E,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(E,[1,2],'FaceColor','r','LineWidth',2);
    plot(E,[1,2],'FaceColor','r','LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

% ------------------------------ END OF CODE ------------------------------
