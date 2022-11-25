function res = testLongDuration_ellipsoid_plot
% testLongDuration_ellipsoid_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check
%
% Syntax:  
%    res = testLongDuration_ellipsoid_plot
%
% Inputs:
%    -
%
% Outputs:
%    res - boolean 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: -

% Author:       Mark Wetzlinger
% Written:      04-August-2020
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

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
    
    % three arguments: object, dimensions, NVpair 'Filled'
    plot(E,[1,2],'Filled',true);
    plot(E,[1,2],'Filled',true,'LineWidth',2);
    plot(E,[1,2],'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(E,[1,2],'r','Filled',true,'LineWidth',2);
    plot(E,[1,2],'r','Filled',true,'LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

if res
    disp('test_plot successful');
else
    disp('test_plot failed');
end

%------------- END OF CODE --------------

