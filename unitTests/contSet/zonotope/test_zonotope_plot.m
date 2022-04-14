function res = test_zonotope_plot
% test_zonotope_plot - unit test function of plot; this function aims
%    to go through many variations of input arguments
%    note: only run-time errors checked, go through manually to check for bugs
%
% Syntax:  
%    res = test_zonotope_plot
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

Z = zonotope(rand(3,5));

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(Z);
    
    % two arguments: object, dimensions
    plot(Z,[1,2]);
    plot(Z,[2,3]);
    
    % three arguments: object, dimensions, linespec
    plot(Z,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(Z,[1,2],'LineWidth',2);
    plot(Z,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Filled'
    plot(Z,[1,2],'Filled',true);
    plot(Z,[1,2],'Filled',true,'LineWidth',2);
    plot(Z,[1,2],'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(Z,[1,2],'r','Filled',true,'LineWidth',2);
    plot(Z,[1,2],'r','Filled',true,'LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
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

