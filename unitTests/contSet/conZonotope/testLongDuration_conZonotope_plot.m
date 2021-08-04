function res = testLongDuration_conZonotope_plot
% testLongDuration_conZonotope_plot - unit test function of plot;
%    this function aims to go through many variations of input arguments
%    note: only run-time errors checked, manual bug check necessary
%
% Syntax:
%    res = testLongDuration_conZonotope_plot
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

% construct zonotope
Z = [0 1.5 -1.5 0.5;0 1 0.5 -1];
A = [1 1 1];
b = 1;
cZono = conZonotope(Z,A,b);

try
    % try all variations in plotting
    figure;
    
    % one argument: object
    plot(cZono);
    
    % two arguments: object, dimensions
    plot(cZono,[1,2]);
    
    % three arguments: object, dimensions, linespec
    plot(cZono,[1,2],'r+');
    
    % three arguments: object, dimensions, NVpairs
    plot(cZono,[1,2],'LineWidth',2);
    plot(cZono,[1,2],'Color',[.6 .6 .6],'LineWidth',2);
    
    % three arguments: object, dimensions, NVpair 'Filled'
    plot(cZono,[1,2],'Filled',true);
    plot(cZono,[1,2],'Filled',true,'LineWidth',2);
    plot(cZono,[1,2],'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % three arguments: object, dimensions, NVpair 'Splits'
    plot(cZono,[1,2],'Splits',4);
    plot(cZono,[1,2],'Splits',4,'LineWidth',2);
    plot(cZono,[1,2],'Splits',4,'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);

    % three arguments: object, dimensions, NVpair 'Template'
    plot(cZono,[1,2],'Template',16);
    plot(cZono,[1,2],'Template',16,'LineWidth',2);
    plot(cZono,[1,2],'Template',16,'Filled',true,'EdgeColor','k','FaceColor',[.8 .8 .8]);
    
    % four arguments: object, dimensions, linespec, NVpairs
    plot(cZono,[1,2],'r','Filled',true,'LineWidth',2);
    plot(cZono,[1,2],'r','Filled',true,'LineWidth',2,'EdgeColor',[.6 .6 .6]);
    
    % close figure
    close;
catch
    close;
    res = false;
end

if res
    disp('testLongDuration_conZonotope_plot successful');
else
    disp('testLongDuration_conZonotope_plot failed');
end

%------------- END OF CODE --------------

